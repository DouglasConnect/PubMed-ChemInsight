import re
import math
import json
import json5
import base64
import datetime
import time
import pandas as pd
import requests
import streamlit as st
from Bio import Entrez
from CompoundResearchHelper import CompoundResearchHelper
from BioInfoRetriever import SynonymRetriever

pd.set_option("display.max_colwidth", 1)

st.set_page_config(page_title="PubMed ChemInsight", page_icon="⚗️")
# Streamlit App Setup 📚
st.markdown(
    """
    <div style="text-align: center;">
        <h1>PubMed ⚗️ ChemInsight \n\n 📑-Unlock Insights from PubMed-📑</h1>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("""
    This tool allows you to search PubMed for Scientific Articles Specific to your Compounds of Interest.
""")


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


# Ensure session state variables exist
if "compounds_text" not in st.session_state:
    st.session_state["compounds_text"] = ""

if "compounds_synonyms_dict" not in st.session_state:
    st.session_state["compounds_synonyms_dict"] = {}


# Convert compounds synonyms dictionary to JSON format
def format_compounds_json(compounds_dict):
    """Convert compound synonyms dictionary into JSON format for display."""
    if compounds_dict:
        return json.dumps(compounds_dict, indent=4, ensure_ascii=False)
    else:
        pass


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

top_syn = st.sidebar.slider(
    "Select the percentage of synonyms to be included in the search", 0, 100, 10
)

# Define layout
col1, col2 = st.columns([3, 1])  # Adjust width ratio as needed
helper = CompoundResearchHelper()

with col1:
    compounds_input = st.text_area(
        "🧪 Enter Compounds (One Chemical Name or CAS Number per Line)",
        value=format_compounds_json(
            st.session_state["compounds_synonyms_dict"]
        ),  # Load JSON as text
        height=150,
        placeholder="Ibuprofen\nMetformin\nAtorvastatin\nOmeprazole",
    )

with col2:
    # Custom CSS for button styling
    st.markdown(
        """
        <style>
            .synonyms-button-container button {
                background-color: #f0f2f6 !important; /* Light Gray */
                color: black !important;
                font-size: 18px !important;
                font-weight: bold !important;
                border-radius: 8px !important;
                cursor: pointer !important;
                transition: all 0.3s !important;
                width: 100% !important;  /* Match text area width */
                height: 80px !important;  /* Adjust height */
                display: flex !important;
                align-items: center !important;
                justify-content: center !important;
                text-align: center !important;
                box-shadow: 2px 4px 6px rgba(0, 0, 0, 0.2) !important;
            }

            .synonyms-button-container button:hover {
                background-color: #c0c0c0 !important; /* Slightly darker Gray */
                box-shadow: 4px 6px 10px rgba(0, 0, 0, 0.3) !important;
            }

            .synonyms-button-container button:active {
                background-color: #a9a9a9 !important; /* Darker Gray */
                box-shadow: 2px 3px 5px rgba(0, 0, 0, 0.4) inset !important;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    # Add a div wrapper to ensure this style applies only to this button
    st.markdown('<div class="synonyms-button-container">', unsafe_allow_html=True)

    # Retrieve Synonyms Button for Compounds
    if st.button(
        "🔍  Retrieve Synonyms for Compounds",
        key="retrieve_compound_synonyms",
        help="Click to Fetch synonyms for compounds",
    ):
        if not compounds_input.strip():
            st.warning("⚠️ Please enter at least one compound!")
        else:
            compounds_list = [
                compound.strip()
                for compound in compounds_input.strip().split("\n")
                if compound.strip()
            ]
            compounds_synonyms_dict = {}
            error_found = False  # Track if an error is found

            for compound in compounds_list:
                compound_name = resolve_compound_name(compound)

                # ✅ Retrieve synonyms and remove empty strings
                synonyms = helper.get_compound_synonyms(compound_name)
                filtered_synonyms = [syn for syn in synonyms if syn.strip()]

                # Calculate the top_syn percentage (rounded up to ensure at least 1 item if the list is small)
                cpd_num_to_keep = max(
                    1, math.ceil(len(filtered_synonyms) * top_syn / 100)
                )
                filtered_synonyms = filtered_synonyms[:cpd_num_to_keep]

                if not filtered_synonyms:
                    st.warning(
                        f"⚠️ No synonyms found for '{compound_name}'. Increase the top synonyms percentage and try again."
                    )

                # ✅ Store the synonyms with the original compound name
                compounds_synonyms_dict[compound_name] = [
                    compound_name
                ] + filtered_synonyms

            # ✅ Store cleaned synonyms dictionary in session state
            st.session_state["compounds_synonyms_dict"] = compounds_synonyms_dict

            # ✅ Immediately update the text area in JSON format
            st.session_state["compounds_text"] = format_compounds_json(
                compounds_synonyms_dict
            )
            st.rerun()  # Force UI refresh

    st.markdown("</div>", unsafe_allow_html=True)  # Close the button container

# # ✅ Display Updated Synonyms in JSON Format
# if "compounds_synonyms_dict" in st.session_state:
#     st.json(st.session_state["compounds_synonyms_dict"])


# Target input with a button beside it
retriever = SynonymRetriever()

# Ensure session state exists
if "targets_text" not in st.session_state:
    st.session_state["targets_text"] = ""

if "synonyms_dict" not in st.session_state:
    st.session_state["synonyms_dict"] = {}


# Convert synonyms_dict to JSON format for text area
def format_synonyms_json(synonyms_dict):
    """Convert dictionary into a pretty JSON format for text area."""
    if synonyms_dict:
        return json.dumps(synonyms_dict, indent=4, ensure_ascii=False)
    else:
        pass


# Define columns
col1, col2 = st.columns([3, 1])  # Adjust width ratio as needed

with col1:
    targets_input = st.text_area(
        "📌 Enter Interaction Targets (One Target per Line in the format of 'target, type')",
        value=format_synonyms_json(
            st.session_state["synonyms_dict"]
        ),  # Load formatted dictionary as text
        height=150,
        placeholder="BRAF, protein\nTP53, gene\naspirin, chemical\nGABA receptor, receptor\nApoptosis, pathway",
    )

with col2:
    # Custom CSS for ONLY this button
    st.markdown(
        """
        <style>
            .element-container:nth-of-type(3) button {
                background-color: #f0f2f6 !important; /* Light Gray */
                color: black !important;
                font-size: 50px !important;
                font-weight: bold !important;
                border-radius: 8px !important;
                border: none !important;
                cursor: pointer !important;
                transition: all 0.3s !important;
                width: 100% !important;  /* Adjust width */
                height: 140px !important;  /* Adjust height */
                display: flex !important;
                align-items: center !important;
                justify-content: center !important;
                text-align: center !important;
            }

            /* Change color on hover */
            .element-container:nth-of-type(3) button:hover {
                background-color: #c0c0c0 !important; /* Slightly darker Gray */
            }

            /* Change color on click */
            .element-container:nth-of-type(3) button:active {
                background-color: #a9a9a9 !important; /* Darker Gray */
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    # Add a div wrapper to ensure this style applies only to this button
    st.markdown('<div class="synonyms-button-container">', unsafe_allow_html=True)

    # Ensure session state exists
    if "error_message" not in st.session_state:
        st.session_state["error_message"] = None
    if "synonyms_dict" not in st.session_state:
        st.session_state["synonyms_dict"] = {}

    # Retrieve Synonyms Button
    if st.button(
        label="🔍  Retrieve Synonyms for\n\n\n Targets",
        key="retrieve_synonyms",
        help="Click to Fetch synonyms for entered targets",
    ):
        if not targets_input.strip():
            st.session_state["error_message"] = "⚠️ Please enter at least one target!"
        else:
            target_list = [
                line.strip().split(",") for line in targets_input.strip().split("\n")
            ]
            synonyms_dict = {}
            error_found = False  # Track if an error is found

            for entry in target_list:
                if len(entry) != 2:  # Ensure exactly two elements exist
                    st.session_state["error_message"] = (
                        f"❌ Invalid format for input: {' '.join(entry)}. "
                        "Use format 'TARGET_NAME, TARGET_TYPE'."
                    )
                    error_found = True
                    continue

                target_name = entry[0].strip()
                target_type = entry[1].strip().lower()

                if target_type not in [
                    "protein",
                    "gene",
                    "chemical",
                    "receptor",
                    "pathway",
                ]:
                    st.session_state["error_message"] = (
                        f"❌ Unknown target type: {target_type}. Choose from "
                        "protein, gene, chemical, receptor, or pathway."
                    )
                    error_found = True
                    continue

                # ✅ Retrieve synonyms inside the loop for each valid target
                synonyms = retriever.get_target_synonyms(target_name, target_type)
                # ✅ Remove empty strings from synonyms list
                filtered_target_synonyms = [syn for syn in synonyms if syn.strip()]
                # Calculate the top_syn percentage (rounded up to ensure at least 1 item if the list is small)
                target_num_to_keep = max(
                    1, math.ceil(len(filtered_target_synonyms) * top_syn / 100)
                )
                filtered_target_synonyms = filtered_target_synonyms[:target_num_to_keep]

                synonyms_dict[target_name] = [
                    target_name
                ] + filtered_target_synonyms  # Add original name

            # ✅ If no errors were found, clear the error message
            if not error_found:
                st.session_state["error_message"] = None

            # ✅ Store synonyms dictionary in session state
            st.session_state["synonyms_dict"] = synonyms_dict
            st.rerun()  # Force Streamlit to refresh and apply updates immediately

    # Display error messages outside the button click block
    if st.session_state["error_message"]:
        st.warning(st.session_state["error_message"])


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


# st.sidebar.write("\n\n\n\n\n\n\n\n\n")
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


def display_summary(compounds, targets):
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

    # compounds_dict = json5.loads(compounds_input)
    # compounds = list(compounds_dict.keys())
    compounds_str = ", ".join(compounds)
    st.markdown(
        f"**Compounds List:** `{compounds_str if compounds_str else 'Not Provided'}`"
    )

    if targets_input:
        # targets_dict = json5.loads(targets_input)
        # targets = list(targets_dict.keys())
        targets_str = ", ".join(targets)
        st.markdown(f"**Interaction Targets List:** `{targets_str}`")
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


def get_key_by_value(dictionary, value):
    for key, val in dictionary.items():
        if val == value:
            return key  # Return the first matching key
    return None  # Return None if not found


def is_valid_json5(text):
    try:
        json5.loads(text)
        return True
    except ValueError:  # json5 raises ValueError instead of JSONDecodeError
        return False


# Main Section for Processing Compounds
if st.button("🚀 Launch Search", help="Click to Start PubMed Search"):
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

    if not compounds_input:
        st.error("Please fill out the compound field.")
    else:
        if is_valid_json5(compounds_input):
            compounds_dict = json5.loads(compounds_input)
            compounds = list(compounds_dict.keys())
        else:
            compounds = [
                compound.strip()
                for compound in compounds_input.split("\n")
                if compound.strip()
            ]
            compounds_dict = {compound: [compound] for compound in compounds}

        # Handle Optional Targets
        if targets_input and targets_input.strip():  # Ensure it's not None before stripping
            if is_valid_json5(targets_input):
                targets_dict = json5.loads(targets_input)
                targets = list(targets_dict.keys())
            else:
                targets = [
                    target.split(",")[0].strip()
                    for target in targets_input.split("\n")
                    if target.strip()
                ]
                targets_dict = {target: [target] for target in targets}
        else:
            targets_dict = {}  # Empty dictionary if no targets provided
            targets = []  # Empty list

        # Display the summary of user-defined configurations
        display_summary(compounds, targets if targets else ["No targets specified"])

        st.info("⏳ Starting article retrieval process...")
        combined_articles = []

        for compound_list in compounds_dict.values():
            if not targets:  # If no targets are specified, search only by compounds
                search_targets = [None]  # Placeholder to allow iteration
            else:
                search_targets = targets_dict.values()

            for target_list in search_targets:
                st.info(
                    f"Searching PubMed for compound: {get_key_by_value(compounds_dict, compound_list)}"
                    + (f" and target: {get_key_by_value(targets_dict, target_list)}" if target_list else "")
                )

                helper.articleList = []
                # Process articles for each search term (compound and optional targets)
                articles_df = helper.process_compound_and_targets(
                    compound_list,
                    target_list if target_list else [],  # Ensure an empty list if no targets
                    start_year,
                    end_year,
                    additional_condition,
                    top_recent_n,
                    bottom_recent_n,
                )

                if not articles_df.empty:
                    st.success(
                        f"✔️ Articles found for: {get_key_by_value(compounds_dict, compound_list)}"
                        + (f" and {get_key_by_value(targets_dict, target_list)}" if target_list else "")
                    )
                    articles_df["url"] = articles_df["url"].str.replace(
                        "https://ncbi.nlm.nih.gov/pubmed/",
                        "https://pubmed.ncbi.nlm.nih.gov/",
                        regex=False,
                    )
                    articles_df["compound"] = get_key_by_value(
                        compounds_dict, compound_list
                    )
                    articles_df["target"] = get_key_by_value(targets_dict, target_list) if target_list else "N/A"
                    articles_df.reset_index(drop=True, inplace=True)
                    st.dataframe(articles_df)
                    combined_articles.append(articles_df)
                else:
                    st.warning(
                        f"No articles found for compound: {get_key_by_value(compounds_dict, compound_list)}"
                        + (f" and {get_key_by_value(targets_dict, target_list)}" if target_list else "")
                    )
                time.sleep(1)

        # Combine articles for all compounds if available
        if len(combined_articles) > 1:
            st.subheader("Combined Articles for All Compounds and Synonyms")
            all_articles_df = pd.concat(
                combined_articles, ignore_index=True
            ).drop_duplicates()
            all_articles_df.reset_index(drop=True, inplace=True)
            st.dataframe(all_articles_df)
            st.info("✅ Done!")
            display_download_button(all_articles_df)
        else:
            st.info("✅ Done!")
            display_download_button(articles_df)
