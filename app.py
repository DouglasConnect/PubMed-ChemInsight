import re
import base64
import datetime
import logging
import os
import time
from io import BytesIO
from itertools import product
from urllib.parse import quote

import pandas as pd
import requests
import streamlit as st
from Bio import Entrez
from metapub import FindIt, PubMedFetcher

pd.set_option("display.max_colwidth", 1)

st.set_page_config(page_title="PubMed ChemInsight", page_icon="🔬")
# Streamlit App Setup 📚
st.markdown(
    """
    <div style="text-align: center;">
        <h1>⚛️ -- PubMed ChemInsight -- ⚛️\nUnlock Insights from PubMed</h1>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("""
    This tool allows you to search PubMed for Scientific Articles Specific to your Compounds of Interest.
""")

# Streamlit Input for User (Main Section)
compounds_input = st.text_area(
    "🧪 Enter Compounds (One Chemical Name or CAS Number per Line)"
)
compounds_list = [
    compound.strip() for compound in compounds_input.split("\n") if compound.strip()
]

genes_input = st.text_area(
    "🧬 Enter Interaction Targets (One Target per Line) (Optional)"
)
genes = [gene.strip() for gene in genes_input.split("\n") if gene.strip()]

# New input for additional keywords from user
additional_keywords_input = st.text_area(
    "🔗 Enter Other Keywords (One Keyword per Line) (Optional)"
)
additional_keywords_list = [
    keyword.strip()
    for keyword in additional_keywords_input.split("\n")
    if keyword.strip()
]

# Modify the additional condition only if additional keywords are provided
additional_condition = (
    f"AND ({' OR '.join([f'{kw}[Title/Abstract]' for kw in additional_keywords_list])})"
    if additional_keywords_list
    else ""
)

file_ = open("images/logo.png", "rb").read()
base64_image = base64.b64encode(file_).decode("utf-8")

st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
st.sidebar.markdown(
    f"""
    <div style="display: flex; align-items: center; justify-content: center;">
        <img src="data:image/png;base64,{base64_image}" alt="Logo" width="200" style="border-radius: 5px;">
    </div>
    """,
    unsafe_allow_html=True,
)

# Sidebar for input fields
st.sidebar.markdown(
    """
    <hr>
    """,
    unsafe_allow_html=True,
)
# st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# st.sidebar.title("Configuration")
# Ask for the user's email
email = st.sidebar.text_input("📧 Enter your email address")
st.sidebar.write("\n\n\n\n\n\n\n\n\n")
top_n = st.sidebar.slider("Select number of top synonyms per compound", 0, 10, 2)
st.sidebar.write("\n\n\n\n\n\n\n\n\n")
top_recent_n = st.sidebar.slider("Select the number of top recent articles", 0, 100, 10)
st.sidebar.write("\n\n\n\n\n\n\n\n\n")
bottom_recent_n = st.sidebar.slider(
    "Select the number of bottom recent articles", 0, 100, 10
)
st.sidebar.write("\n\n\n\n\n\n\n\n\n")
start_year = st.sidebar.number_input(
    "Start Year",
    value=2000,
    step=1,
    min_value=1900,
    max_value=datetime.datetime.now().year,
)
st.sidebar.write("\n\n\n\n\n\n\n\n\n")
end_year = st.sidebar.number_input(
    "End Year",
    value=datetime.datetime.now().year,
    step=1,
    min_value=1900,
    max_value=datetime.datetime.now().year,
)

# Sidebar for Configuration
st.sidebar.markdown(
    """
    <hr>
    <div style="text-align: center;">
        <p style="font-size: 12px; color: gray;">
            © 2024 Edelweiss Connect - Developed by Asmaa A. Abdelwahab
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.sidebar.markdown(
    """
    <div style="display: flex; align-items: center; justify-content: center;">
        <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="text-decoration: none;">
            <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" alt="GitHub Logo" style="width:40px; height:40px; margin-right: 10px;">
        </a>
        <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="text-decoration: none;">
            <p style="font-size: 16px; font-weight: bold; color: black; margin: 0;">@asmaa-a-abdelwahab</p>
        </a>
    </div>
    """,
    unsafe_allow_html=True,
)


# Cache to prevent redundant queries
@st.cache_data(show_spinner=False)
def cache_synonyms(compound):
    return CompoundResearchHelper().most_common_synonyms(compound, top_n)


def cas_to_iupac(cas_number):
    """
    Converts a CAS number to an IUPAC name using the CACTUS Chemical Identifier Resolver service.

    Parameters:
    cas_number (str): The CAS number to be converted.

    Returns:
    str: The corresponding IUPAC name or an error message if the conversion fails.
    """
    try:
        # Construct the URL for the CACTUS service with CAS number and request for IUPAC name
        url = f"https://cactus.nci.nih.gov/chemical/structure/{cas_number}/iupac_name"

        # Send a GET request to the URL
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            return response.text.strip()
        else:
            return f"Error: Unable to retrieve IUPAC name (HTTP Status Code: {response.status_code})"
    except Exception as e:
        return f"An error occurred: {str(e)}"


def is_cas_number(compound):
    """Check if the provided string matches the CAS number format."""
    return bool(re.match(r"^\d{2,7}-\d{2}-\d$", compound))


def resolve_compound_name(compound):
    """
    Resolve a compound name. If the compound is a CAS number, convert it to the IUPAC name.

    Parameters:
    compound (str): The compound input by the user.

    Returns:
    str: The resolved compound name.
    """
    if is_cas_number(compound):
        iupac_name = cas_to_iupac(compound)
        if "Error" in iupac_name or "An error occurred" in iupac_name:
            st.warning(f"Failed to convert CAS number '{compound}': {iupac_name}")
            return compound  # Return the original CAS number if conversion fails
        else:
            st.info(f"CAS number '{compound}' converted to IUPAC name '{iupac_name}'")
            return iupac_name
    else:
        return compound


class CompoundResearchHelper:
    """A class to fetch compound synonyms and retrieve relevant articles from PubMed."""

    def __init__(self, retmax=1000):
        self.pubmed = PubMedFetcher()
        self.retmax = retmax
        self.synonym_data = {}
        self.articleList = []
        self.current_year = datetime.datetime.now().year

    def get_compound_synonyms(self, compound_name):
        """Retrieve compound synonyms from PubChem."""
        try:
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
            cid_response = requests.get(cid_url)
            cid_response.raise_for_status()

            cid_data = cid_response.json()
            cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]
            if not cid:
                logging.warning(f"No CID found for {compound_name}")
                return []

            synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url)
            synonyms_response.raise_for_status()

            synonyms_data = synonyms_response.json()
            all_synonyms = (
                synonyms_data.get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )
            return all_synonyms or []

        except requests.HTTPError as http_err:
            logging.error(f"HTTP error occurred: {http_err}")
        except Exception as err:
            logging.error(f"An error occurred: {err}")
        return []

    def get_synonym_frequency(self, synonyms):
        """Check the frequency of each synonym in PubMed literature."""
        synonym_frequency = {}
        for synonym in synonyms:
            try:
                encoded_synonym = quote(f'"{synonym}"')
                handle = Entrez.esearch(
                    db="pubmed",
                    term=encoded_synonym,
                    retmax=0,
                )
                record = Entrez.read(handle)
                synonym_frequency[synonym] = int(record.get("Count", 0))
                handle.close()
                time.sleep(0.3)  # Respect PubMed rate limit
            except Exception as e:
                logging.error(f"Error fetching frequency for {synonym}: {e}")
                continue
        return synonym_frequency

    def most_common_synonyms(self, compound_name, top_n=5):
        """Get the most commonly used synonyms based on literature frequency."""
        if compound_name in self.synonym_data:
            return [synonym for synonym, _ in self.synonym_data[compound_name][:top_n]]

        all_synonyms = self.get_compound_synonyms(compound_name)
        if not all_synonyms:
            return []

        frequency = self.get_synonym_frequency(all_synonyms)
        sorted_synonyms = sorted(
            frequency.items(), key=lambda item: item[1], reverse=True
        )
        self.synonym_data[compound_name] = sorted_synonyms
        return [compound_name] + [synonym for synonym, _ in sorted_synonyms[:top_n]]

    def fetch_articles(self, search_term, retmax=1000, start_year=2000, end_year=None):
        """Fetch articles' metadata from PubMed based on a given search term and date range."""
        date_range = (
            f' AND ("{start_year}/01/01"[PDat] : "{end_year}/12/31"[PDat])'
            if end_year
            else ""
        )
        full_search_term = f"{search_term}{date_range}"

        try:
            pmids = self.pubmed.pmids_for_query(
                full_search_term, retmax=retmax, sort="relevance"
            )
        except Exception as e:
            logging.error(f"Error fetching pmids for {search_term}: {e}")
            return []

        articles = []
        for pmid in pmids:
            try:
                article = self.pubmed.article_by_pmid(pmid)
                article_dict = article.to_dict()
                keys_to_keep = [
                    "title",
                    "pmid",
                    "url",
                    "authors",
                    "doi",
                    "pmc",
                    "issn",
                    "mesh",
                    "chemicals",
                    "journal",
                    "abstract",
                    "year",
                ]
                article_dict = {k: article_dict.get(k, None) for k in keys_to_keep}
                articles.append(article_dict)
            except Exception as e:
                logging.error(f"Error processing article with pmid {pmid}: {e}")
        return articles

    def filter_articles_by_recent(self, df, top_n, bottom_n):
        """Filter the top and bottom recent articles based on user input."""
        if "year" in df.columns:
            df["year"] = pd.to_numeric(df["year"], errors="coerce")

        # Sort the dataframe by the 'year' column in descending order
        df_sorted = df.sort_values(by="year", ascending=False)

        # Convert unhashable columns (like lists and dictionaries) to strings
        df_sorted = df_sorted.applymap(
            lambda x: str(x) if isinstance(x, (list, dict)) else x
        )

        # Check for zero values and filter accordingly
        top_recent = df_sorted.head(top_n) if top_n > 0 else pd.DataFrame()
        bottom_recent = df_sorted.tail(bottom_n) if bottom_n > 0 else pd.DataFrame()

        # Concatenate the results and drop duplicates based on specific columns (e.g., title and pmid)
        filtered_df = pd.concat([top_recent, bottom_recent]).drop_duplicates(
            subset=["title", "pmid"]
        )

        return filtered_df

    def process_compound_and_genes(
        self,
        compound,
        genes,
        start_year,
        end_year,
        additional_condition,
        top_n,
        bottom_n,
    ):
        """Process a single compound and multiple interaction targets for synonym lookup and article fetching."""
        logging.info(f"Processing compound: {compound}")

        # Check if synonyms are enabled by checking if top_n is greater than zero
        if top_n > 0:
            top_synonyms = cache_synonyms(compound)
        else:
            # If top_n is zero, use only the original compound name
            top_synonyms = [compound]

        # Handle queries with or without interaction targets (genes)
        if genes:
            queries = [
                f"({synonym}[Title/Abstract]) AND ({gene}[Title/Abstract]) {additional_condition}"
                for synonym, gene in product(top_synonyms, genes)
            ]
        else:
            queries = [
                f"({synonym}[Title/Abstract]) {additional_condition}"
                for synonym in top_synonyms
            ]

        # Fetch articles for each combination of synonym and interaction target
        for query in queries:
            self.articleList.extend(
                self.fetch_articles(
                    query, retmax=self.retmax, start_year=start_year, end_year=end_year
                )
            )

        # Create a DataFrame from the articles and filter them based on the user's input
        if self.articleList:
            df = pd.DataFrame(self.articleList)
            df_filtered = self.filter_articles_by_recent(df, top_n, bottom_n)
            return df_filtered
        return pd.DataFrame()


# Run when the user clicks the "Search" button
def cas_to_iupac_pubchem(cas_number, retries=3, backoff_factor=2):
    """
    Converts a CAS number to an IUPAC name using the PubChem PUG-REST API.
    Implements a retry mechanism with exponential backoff for handling temporary server errors (503).

    Parameters:
    cas_number (str): The CAS number to be converted.
    retries (int): Number of times to retry the request in case of server errors.
    backoff_factor (int): Factor by which the wait time increases after each retry.

    Returns:
    str: The corresponding IUPAC name or an error message if the conversion fails.
    """
    try:
        # Construct the URL for the PubChem PUG-REST API to retrieve CID
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/cids/JSON"

        # Retry loop
        for attempt in range(retries):
            try:
                cid_response = requests.get(cid_url)
                cid_response.raise_for_status()  # Raise HTTPError for bad responses (4xx and 5xx)
                cid_data = cid_response.json()
                cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]

                if not cid:
                    return f"Error: No CID found for CAS number '{cas_number}'"

                # Fetch the IUPAC name using the CID
                iupac_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                iupac_response = requests.get(iupac_url)
                iupac_response.raise_for_status()
                iupac_data = iupac_response.json()
                iupac_name = (
                    iupac_data.get("PropertyTable", {})
                    .get("Properties", [{}])[0]
                    .get("IUPACName", None)
                )

                if iupac_name:
                    return iupac_name
                else:
                    return "Error: No IUPAC name found"

            except requests.exceptions.HTTPError as http_err:
                if cid_response.status_code == 503:  # Server busy or unavailable
                    wait_time = backoff_factor**attempt
                    st.warning(
                        f"PubChem service temporarily unavailable. Retrying in {wait_time} seconds..."
                    )
                    time.sleep(wait_time)  # Wait before retrying
                else:
                    return f"HTTP error occurred: {http_err}"

        # If all retries fail, return a final error message
        return "Error: PubChem service unavailable after multiple attempts."

    except Exception as err:
        return f"An error occurred: {str(err)}"


def resolve_compound_name(compound):
    """
    Resolve a compound name. If the compound is a CAS number, attempt to convert it to the IUPAC name.
    It tries PubChem first and falls back to CACTUS if needed.

    Parameters:
    compound (str): The compound input by the user.

    Returns:
    str: The resolved compound name.
    """
    if is_cas_number(compound):
        # Try converting using PubChem first
        iupac_name = cas_to_iupac_pubchem(compound)
        if "Error" in iupac_name or "An error occurred" in iupac_name:
            st.warning(
                f"PubChem failed to convert CAS number '{compound}': {iupac_name}"
            )
            # Try converting using CACTUS as a fallback
            iupac_name = cas_to_iupac(compound)
            if "Error" in iupac_name or "An error occurred" in iupac_name:
                st.warning(
                    f"CACTUS also failed to convert CAS number '{compound}': {iupac_name}"
                )
                return compound  # Return the original CAS number if conversion fails
            else:
                st.info(
                    f"CAS number '{compound}' converted to IUPAC name '{iupac_name}' using CACTUS"
                )
                return iupac_name
        else:
            st.info(
                f"CAS number '{compound}' converted to IUPAC name '{iupac_name}' using PubChem"
            )
            return iupac_name
    else:
        return compound


def display_summary():
    """
    Display a summary of the user-defined configurations.
    """
    st.markdown("### Summary of Your Selections")

    st.markdown(f"**Email Address:** `{email if email else 'Not Provided'}`")

    compounds_str = ", ".join(compounds_list)
    st.markdown(
        f"**Compounds List:** `{compounds_str if compounds_str else 'Not Provided'}`"
    )

    if genes:
        genes_str = ", ".join(genes)
        st.markdown(f"**Interaction Targets List:** `{genes_str}`")
    else:
        st.markdown("**Interaction Targets List:** `Not Provided`")

    if additional_keywords_list:
        keywords_str = ", ".join(additional_keywords_list)
        st.markdown(f"**Additional Keywords:** `{keywords_str}`")
    else:
        st.markdown("**Additional Keywords:** `Not Provided`")

    # Show the number of synonyms to be retrieved or indicate that no synonyms are selected
    if top_n > 0:
        st.markdown(f"**Top Synonyms per Compound:** `{top_n}`")
    else:
        st.markdown(f"**Top Synonyms per Compound:** `Not Selected`")

    st.markdown(
        f"**Number of Top Recent Articles:** `{top_recent_n if top_recent_n > 0 else 'Not Selected'}`"
    )
    st.markdown(
        f"**Number of Bottom Recent Articles:** `{bottom_recent_n if bottom_recent_n > 0 else 'Not Selected'}`"
    )
    st.markdown(f"**Year Range:** `{start_year} to {end_year}`")

    st.markdown("---")  # A horizontal line to separate the summary from the results


# Main Section for Processing Compounds
if st.button("🚀 Launch Search"):
    # Validate the email input using basic checks
    if email:
        if "@" not in email or "." not in email.split("@")[-1]:
            st.error("Please enter a valid email address.")
        else:
            st.success(f"Email address '{email}' is valid!")

    # Ensure the email is set for Entrez
    Entrez.email = email

    if not compounds_list:
        st.error("Please fill out the compound field.")
    else:
        # Display the summary of user-defined configurations
        display_summary()

        st.info("⏳ Starting article retrieval process...")
        helper = CompoundResearchHelper()

        all_articles = []
        for compound in compounds_list:
            resolved_compound = resolve_compound_name(compound)
            st.info(f"Processing compound: {resolved_compound}")

            helper.articleList = []

            articles_df = helper.process_compound_and_genes(
                resolved_compound,
                genes,
                start_year,
                end_year,
                additional_condition,
                top_recent_n,
                bottom_recent_n,
            )

            if not articles_df.empty:
                articles_df.reset_index(drop=True, inplace=True)

                output = BytesIO()
                articles_df.to_csv(output, index=False)
                output.seek(0)

                filename = f"{resolved_compound}_filtered_articles.csv"
                with open(filename, "wb") as f:
                    f.write(output.read())

                st.success(f"✔️ Articles processed for compound: {resolved_compound}")
                st.dataframe(articles_df)

                all_articles.append(articles_df)
            else:
                st.warning(f"No articles found for compound: {resolved_compound}")

        if all_articles:
            combined_df = pd.concat(all_articles, ignore_index=True)
            combined_df.drop_duplicates(
                subset=["title", "pmid"], keep="first", inplace=True
            )
            combined_df.to_csv("combined_pubmed_articles.csv", index=False)

            st.subheader("Combined Articles for All Compounds")
            st.dataframe(combined_df)
            st.success(f"✔️ Combined articles saved to 'combined_pubmed_articles.csv'")
        else:
            st.warning("No articles found for any of the compounds.")
        st.info("✅ Done!")
