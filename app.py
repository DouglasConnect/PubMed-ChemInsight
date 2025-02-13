import re
import base64
import datetime
import time
import pandas as pd
import requests
import streamlit as st
from Bio import Entrez
from CompoundResearchHelper import CompoundResearchHelper

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

st.sidebar.markdown(
    f"""
    <div style="display: flex; align-items: center; justify-content: center; padding-bottom: 10px;">
        <!-- Logo -->
        <div style="display: flex; align-items: center; margin-right: 10px;">
            <img src="data:image/png;base64,{base64_image}" alt="Logo" width="120" style="border-radius: 5px;">
        </div>
        <!-- Separator -->
        <div style="width: 4px; height: 30px; background-color: #ccc; margin-right: 10px;"></div>
        <!-- Text -->
        <div style="text-align: center;">
            <a href="https://doi.org/10.5281/zenodo.14771565">
                <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14771565.svg" alt="DOI">
            </a>
        </div>
        <hr>
    </div>
    """,
    unsafe_allow_html=True,
)
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# # Sidebar for logo
# st.sidebar.markdown(
#     f"""
#     <div style="display: flex; align-items: center; justify-content: center;">
#         <img src="data:image/png;base64,{base64_image}" alt="Logo" width="150" style="border-radius: 5px;">
#     </div>
#     """,
#     unsafe_allow_html=True,
# )

# # Sidebar for input fields
# st.sidebar.markdown(
#     """
#     <div style="text-align: center;">
#         <a href="https://doi.org/10.5281/zenodo.14771565">
#             <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14771565.svg" alt="DOI">
#         </a>
#     </div>
#     <hr>
#     """,
#     unsafe_allow_html=True,
# )
# st.sidebar.markdown(
#     """
#     [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14774229.svg)](https://doi.org/10.5281/zenodo.14774229)
#     """,
#     unsafe_allow_html=True,
# )
# st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# st.sidebar.title("Configuration")
# Ask for the user's email
email = st.sidebar.text_input("📧 Enter your email address (Preferred)")
api_key = st.sidebar.text_input(
    "🔑 Enter your NCBI API key (Preferred)", type="password"
)
st.sidebar.markdown(
    """
    <div style="text-align: center;">
        <p style="font-size:12px; color:black;">
            <a href="https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us" target="_blank" style="font-size:14px; color:black;">Create NCBI API key for a better performance</a>
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)
# st.sidebar.write("\n\n\n\n\n\n\n\n\n")
top_syn = st.sidebar.slider("Select number of top synonyms per compound", 0, 10, 2)
top_recent_n = st.sidebar.slider("Select the number of top recent articles", 0, 100, 10)
bottom_recent_n = st.sidebar.slider(
    "Select the number of bottom recent articles", 0, 100, 10
)
# st.sidebar.write("\n\n\n\n\n\n\n\n\n")
start_year = st.sidebar.number_input(
    "Start Year",
    value=2000,
    step=1,
    min_value=1900,
    max_value=datetime.datetime.now().year,
)
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
            <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" alt="GitHub Logo" style="width:30px; height:30px; margin-right: 10px;">
        </a>
        <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="text-decoration: none;">
            <p style="font-size: 14px; font-weight: bold; color: black; margin: 0;">@asmaa-a-abdelwahab</p>
        </a>
    </div>
    """,
    unsafe_allow_html=True,
)


def is_cas_number(compound):
    """Check if the provided string matches the CAS number format."""
    return bool(re.match(r"^\d{2,7}-\d{2}-\d$", compound))


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

        time.sleep(0.2)
        # If all retries fail, return a final error message
        return "Error: PubChem service unavailable after multiple attempts."

    except Exception as err:
        return f"An error occurred: {str(err)}"


def cas_to_iupac(cas_number):
    """
    Converts a CAS number to an IUPAC name using the CACTUS server.

    Parameters:
    cas_number (str): The CAS number to be converted.

    Returns:
    str: The corresponding IUPAC name or an error message if the conversion fails.
    """

    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{cas_number}/iupac_name"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
        else:
            return f"Error: Unable to retrieve IUPAC name (HTTP Status Code: {response.status_code})"
    except Exception as e:
        return f"An error occurred: {str(e)}"


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
    Displays a summary of the user's selections including email, API key, compounds, interaction targets,
    additional keywords, number of recent articles, and year range.

    Utilizes Streamlit to present the information in a markdown format.

    Parameters:
    None

    Returns:
    None
    """

    st.markdown("### Summary of Your Selections")

    st.markdown(f"**Email Address:** `{email if email else 'Not Provided'}`")
    st.markdown(f"**NCBI API Key:** `{api_key if api_key else 'Not Provided'}`")

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

    st.markdown(
        f"**Number of Top Recent Articles:** `{top_recent_n if top_recent_n > 0 else 'Not Selected'}`"
    )
    st.markdown(
        f"**Number of Bottom Recent Articles:** `{bottom_recent_n if bottom_recent_n > 0 else 'Not Selected'}`"
    )
    st.markdown(f"**Year Range:** `{start_year} to {end_year}`")

    st.markdown("---")  # A horizontal line to separate the summary from the results


@st.fragment
def display_download_button(all_articles_df):
    # Convert DataFrame to CSV in-memory
    csv = all_articles_df.to_csv(index=False).encode("utf-8")

    # Add a download button
    if st.download_button(
        label="📥 Download Combined Articles CSV",
        data=csv,
        file_name="combined_pubmed_articles.csv",
        mime="text/csv",
    ):
        st.success("✔️ Combined articles saved to 'combined_pubmed_articles.csv'")


# Main Section for Processing Compounds
if st.button("🚀 Launch Search"):
    # Validate the email input using basic checks
    if email:
        if "@" not in email or "." not in email.split("@")[-1]:
            st.error("Please enter a valid email address.")
        else:
            st.success(f"Email address '{email}' is valid!")

    # Set the email for Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    if not compounds_list:
        st.error("Please fill out the compound field.")
    else:
        # Display the summary of user-defined configurations
        display_summary()

        st.info("⏳ Starting article retrieval process...")
        helper = CompoundResearchHelper()

        combined_articles = []
        for compound in compounds_list:
            st.info(f"Processing compound: {compound}")
            resolved_compound = resolve_compound_name(compound)
            st.info(f"Searching PubMed for compound: {resolved_compound}")

            helper.articleList = []
            # Process articles for each search term (compound and synonyms)
            articles_df = helper.process_compound_and_genes(
                resolved_compound,
                genes,
                start_year,
                end_year,
                additional_condition,
                top_recent_n,
                bottom_recent_n,
                top_syn,
            )

            if not articles_df.empty:
                st.success(f"✔️ Articles found for: {compound}")
                articles_df["compound"] = compound
                articles_df.reset_index(
                    drop=True, inplace=True
                )  # Tag articles with compound name
                st.dataframe(articles_df)
                combined_articles.append(articles_df)
            else:
                st.warning(f"No articles found for compound: {compound}")
            time.sleep(1)

        # Combine articles for all compounds if available
        if combined_articles:
            st.subheader("Combined Articles for All Compounds and Synonyms")
            all_articles_df = pd.concat(
                combined_articles, ignore_index=True
            ).drop_duplicates()
            all_articles_df.reset_index(drop=True, inplace=True)
            st.dataframe(all_articles_df)
            st.info("✅ Done!")
            display_download_button(all_articles_df)
        else:
            st.warning("No articles found for any of the compounds.")
