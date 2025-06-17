import re
import ast
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
import os
import logging
import threading
import queue
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

# Set pandas display option
pd.set_option("display.max_colwidth", 1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("app_debug.log", mode="a"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# Email configuration (set these in your environment variables)
EMAIL_ADDRESS = "asmaa@edelweissconnect.com"
EMAIL_PASSWORD = "tkkn pawl qzrb fbke"

# Streamlit App Setup
st.set_page_config(page_title="PubMed ChemInsight", page_icon="⚗️")
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


# Function to send email with optional attachment
def send_email(to_email, subject, body, attachment_path=None):
    try:
        msg = MIMEMultipart()
        msg["From"] = EMAIL_ADDRESS
        msg["To"] = to_email
        msg["Subject"] = subject
        msg.attach(MIMEText(body, "plain"))

        if attachment_path and os.path.exists(attachment_path):
            with open(attachment_path, "rb") as attachment:
                part = MIMEBase("application", "octet-stream")
                part.set_payload(attachment.read())
                encoders.encode_base64(part)
                part.add_header(
                    "Content-Disposition",
                    f"attachment; filename={os.path.basename(attachment_path)}",
                )
                msg.attach(part)

        with smtplib.SMTP("smtp.gmail.com", 587) as server:
            server.starttls()
            server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
            server.send_message(msg)
        logger.info(f"Email sent to {to_email}")
    except Exception as e:
        logger.error(f"Error sending email: {e}")


def safe_parse_publication_types(x):
    try:
        if isinstance(x, dict):
            # Handle nested dictionaries or unexpected formats
            values = []
            for v in x.values():
                if isinstance(v, (list, dict)):
                    values.append(str(v))
                else:
                    values.append(v)
            return ", ".join(values)
        elif isinstance(x, str) and x.strip().startswith("{"):
            parsed = ast.literal_eval(x)
            if isinstance(parsed, dict):
                values = []
                for v in parsed.values():
                    if isinstance(v, (list, dict)):
                        values.append(str(v))
                    else:
                        values.append(v)
                return ", ".join(values)
            else:
                return x
        else:
            return x
    except Exception as e:
        logging.warning(f"⚠️ Failed to parse publication_types: {x}, error: {e}")
        return str(x)  # Convert to string as a fallback


# Function to perform PubMed search in the background
def perform_pubmed_search(task):
    try:
        Entrez.email = task["email"]
        if task["api_key"]:
            Entrez.api_key = task["api_key"]

        helper = CompoundResearchHelper()
        additional_condition = (
            f"AND ({' OR '.join([f'{kw}[Title/Abstract]' for kw in task['additional_keywords_list']])})"
            if task["additional_keywords_list"]
            else ""
        )

        combined_articles = []

        for compound_original, compound_synonyms in task["compounds_dict"].items():
            if not task["targets_dict"]:
                target_pairs = [(None, None)]
            else:
                target_pairs = list(task["targets_dict"].items())

            for target_original, target_synonyms in target_pairs:
                logger.info(
                    f"🔎 Searching PubMed for compound: {compound_original}"
                    + (f" and target: {target_original}" if target_original else "")
                )

                helper.articleList = []
                articles_df = helper.process_compound_and_targets(
                    compounds=compound_synonyms,
                    genes=target_synonyms if target_synonyms else [],
                    start_year=task["start_year"],
                    end_year=task["end_year"],
                    additional_condition=additional_condition,
                    n_articles=task["n_articles_per_pair"],
                    article_type_query=task["article_type_query"],
                )

                if not articles_df.empty:
                    # Fix URL formatting
                    articles_df["url"] = articles_df["url"].str.replace(
                        "https://ncbi.nlm.nih.gov/pubmed/",
                        "https://pubmed.ncbi.nlm.nih.gov/",
                        regex=False,
                    )

                    articles_df["compound"] = get_key_by_value(
                        task["compounds_dict"], compound_synonyms
                    )
                    articles_df["target"] = (
                        get_key_by_value(task["targets_dict"], target_synonyms)
                        if target_synonyms
                        else "N/A"
                    )
                    articles_df.reset_index(drop=True, inplace=True)
                    articles_df["publication_types"] = articles_df[
                        "publication_types"
                    ].apply(safe_parse_publication_types)

                    # Convert unhashable types (lists and dictionaries) to strings
                    for col in articles_df.columns:
                        # Check for lists
                        if articles_df[col].apply(lambda x: isinstance(x, list)).any():
                            articles_df[col] = articles_df[col].apply(
                                lambda x: ", ".join(map(str, x))
                                if isinstance(x, list)
                                else x
                            )
                        # Check for dictionaries
                        if articles_df[col].apply(lambda x: isinstance(x, dict)).any():
                            articles_df[col] = articles_df[col].apply(
                                lambda x: ", ".join(f"{k}: {v}" for k, v in x.items())
                                if isinstance(x, dict)
                                else x
                            )

                    combined_articles.append(articles_df)
                else:
                    logger.warning(
                        f"⚠️ No articles found for compound: {compound_original}"
                        + (f" and target: {target_original}" if target_original else "")
                    )

        # Determine if synonyms were retrieved
        compounds_have_synonyms = any(
            len(synonyms) > 1 for synonyms in task["compounds_dict"].values()
        )
        targets_have_synonyms = any(
            len(synonyms) > 1 for synonyms in task["targets_dict"].values()
        )

        # Generate synonym retrieval message
        if compounds_have_synonyms and targets_have_synonyms:
            synonym_message = "Synonyms were retrieved for both compounds and targets."
        elif compounds_have_synonyms:
            synonym_message = "Synonyms were retrieved for compounds only."
        elif targets_have_synonyms:
            synonym_message = "Synonyms were retrieved for targets only."
        else:
            synonym_message = "No synonyms were retrieved for compounds or targets."

        # Extract compounds and targets for the summary
        compounds = list(task["compounds_dict"].keys())
        targets = (
            list(task["targets_dict"].keys())
            if task["targets_dict"]
            else ["No targets specified"]
        )

        # Generate the summary
        summary = generate_summary(task, compounds, targets)

        # Combine the email body
        if combined_articles:
            all_articles_df = pd.concat(
                combined_articles, ignore_index=True
            ).drop_duplicates()
            all_articles_df.reset_index(drop=True, inplace=True)

            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            csv_path = f"pubmed_chminsight_results_{timestamp}.csv"

            all_articles_df.to_csv(csv_path, index=False)

            # Create the email body with summary and synonym info
            email_body = (
                "Thank you for using PubMed ChemInsight!\n\n"
                "Please find attached the CSV file with your PubMed search results.\n\n"
                f"{synonym_message}\n\n"
                f"{summary}"
            )

            send_email(
                to_email=task["email"],
                subject="Your PubMed Search Results",
                body=email_body,
                attachment_path=csv_path,
            )
            if os.path.exists(csv_path):
                os.remove(csv_path)
        else:
            # Create the email body for the "no articles found" case
            email_body = (
                "Thank you for using PubMed ChemInsight!\n\n"
                "No articles were found for your search criteria.\n\n"
                f"{synonym_message}\n\n"
                f"{summary}"
            )

            send_email(
                to_email=task["email"],
                subject="PubMed Search Results",
                body=email_body,
            )

    except Exception as e:
        logger.error(f"❌ Error during PubMed search: {e}")


# Background worker function
def worker(task_queue):
    while True:
        task = task_queue.get()
        if task is None:  # Exit condition
            break
        try:
            perform_pubmed_search(task)
        except Exception as e:
            logger.error(f"Error processing task: {e}")
        finally:
            task_queue.task_done()


# Initialize task queue and worker thread
@st.cache_resource
def get_task_queue():
    task_queue = queue.Queue()
    worker_thread = threading.Thread(target=worker, args=(task_queue,), daemon=True)
    worker_thread.start()
    return task_queue


task_queue = get_task_queue()


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


# Sidebar setup (unchanged except for logo handling)
file_ = open("images/logo.png", "rb").read()
base64_image = base64.b64encode(file_).decode("utf-8")
st.sidebar.markdown(
    f"""
    <div style="display: flex; align-items: center; justify-content: center; padding-bottom: 10px;">
        <img src="data:image/png;base64,{base64_image}" alt="Logo" width="120" style="border-radius: 5px;">
        <div style="width: 4px; height: 30px; background-color: #ccc; margin-right: 10px;"></div>
        <a href="https://doi.org/10.5281/zenodo.14771565">
            <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14771565.svg" alt="DOI">
        </a>
        <hr>
    </div>
    """,
    unsafe_allow_html=True,
)

# Sidebar inputs (unchanged)
email = st.sidebar.text_input("📧 Enter your email address (Required)")
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

# New inputs for synonyms and articles
num_synonyms_per_compound = st.sidebar.number_input(
    "Number of Synonyms per Compound", min_value=0, value=5, step=1
)
num_synonyms_per_target = st.sidebar.number_input(
    "Number of Synonyms per Target", min_value=0, value=5, step=1
)
n_articles_per_pair = st.sidebar.number_input(
    "Number of Articles per Compound-Target Pair", min_value=1, value=10, step=1
)

# Session state initialization (unchanged)
if "compounds_synonyms_dict" not in st.session_state:
    st.session_state["compounds_synonyms_dict"] = {}
if "resolved_compounds" not in st.session_state:
    st.session_state["resolved_compounds"] = {}


def format_compounds_json(compounds_dict):
    return (
        json.dumps(compounds_dict, indent=4, ensure_ascii=False)
        if compounds_dict
        else ""
    )


# Compound input section (unchanged)
col1, col2 = st.columns([3, 1])
helper = CompoundResearchHelper()

with col1:
    compounds_input = st.text_area(
        "🧪 Enter Compounds (One Chemical Name or CAS Number per Line)",
        height=150,
        placeholder="Ibuprofen\nMetformin\nAtorvastatin\nOmeprazole",
        value=st.session_state.get("formatted_compounds", ""),
        key="compounds_input",
    )
    if not st.session_state["compounds_input"].strip():
        st.session_state["compounds_synonyms_dict"] = {}
        st.session_state["resolved_compounds"] = {}

with col2:
    st.markdown(
        """
        <style>
            .synonyms-button-container button {
                background-color: #f0f2f6 !important;
                color: black !important;
                font-size: 18px !important;
                font-weight: bold !important;
                border-radius: 8px !important;
                cursor: pointer !important;
                transition: all 0.3s !important;
                width: 100% !important;
                height: 80px !important;
                display: flex !important;
                align-items: center !important;
                justify-content: center !important;
                text-align: center !important;
                box-shadow: 2px 4px 6px rgba(0, 0, 0, 0.2) !important;
            }
            .synonyms-button-container button:hover {
                background-color: #c0c0c0 !important;
                box-shadow: 4px 6px 10px rgba(0, 0, 0, 0.3) !important;
            }
            .synonyms-button-container button:active {
                background-color: #a9a9a9 !important;
                box-shadow: 2px 3px 5px rgba(0, 0, 0, 0.4) inset !important;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown('<div class="synonyms-button-container">', unsafe_allow_html=True)
    if st.button(
        "🔍 Retrieve Synonyms for Compounds", key="retrieve_compound_synonyms"
    ):
        compounds_list = [
            compound.strip()
            for compound in (st.session_state["compounds_input"] or "")
            .strip()
            .split("\n")
            if compound.strip()
        ]
        if not compounds_list:
            st.warning("⚠️ Please enter at least one compound!")
        else:
            resolved_compounds = {
                compound: resolve_compound_name(compound) for compound in compounds_list
            }
            st.session_state["resolved_compounds"] = resolved_compounds
            compounds_synonyms_dict = {}
            for original_name, resolved_name in resolved_compounds.items():
                try:
                    synonyms = helper.get_compound_synonyms(resolved_name)
                    filtered_synonyms = [syn for syn in synonyms if syn.strip()]
                    filtered_synonyms = filtered_synonyms[:num_synonyms_per_compound]
                    if not filtered_synonyms:
                        st.warning(
                            f"⚠️ No synonyms found for '{resolved_name}'. Only the original name will be used."
                        )
                        compounds_synonyms_dict[original_name] = [resolved_name]
                    else:
                        compounds_synonyms_dict[original_name] = [
                            resolved_name
                        ] + filtered_synonyms
                except Exception as e:
                    st.error(f"Error retrieving synonyms for '{resolved_name}': {e}")
            st.session_state["compounds_synonyms_dict"] = compounds_synonyms_dict
            st.session_state["formatted_compounds"] = format_compounds_json(
                compounds_synonyms_dict
            )
            st.rerun()
    st.markdown("</div>", unsafe_allow_html=True)

# Target input section (unchanged)
retriever = SynonymRetriever()
if "targets_text" not in st.session_state:
    st.session_state["targets_text"] = ""
if "synonyms_dict" not in st.session_state:
    st.session_state["synonyms_dict"] = {}


def format_synonyms_json(synonyms_dict):
    return (
        json.dumps(synonyms_dict, indent=4, ensure_ascii=False) if synonyms_dict else ""
    )


col1, col2 = st.columns([3, 1])
with col1:
    targets_input = st.text_area(
        "📌 Enter Interaction Targets (One Target per Line in the format of 'target, type')",
        value=format_synonyms_json(st.session_state["synonyms_dict"]),
        height=150,
        placeholder="BRAF, protein\nTP53, gene\naspirin, chemical\nGABA receptor, receptor\nApoptosis, pathway",
    )
with col2:
    st.markdown(
        """
        <style>
            .element-container:nth-of-type(3) button {
                background-color: #f0f2f6 !important;
                color: black !important;
                font-size: 50px !important;
                font-weight: bold !important;
                border-radius: 8px !important;
                border: none !important;
                cursor: pointer !important;
                transition: all 0.3s !important;
                width: 100% !important;
                height: 140px !important;
                display: flex !important;
                align-items: center !important;
                justify-content: center !important;
                text-align: center !important;
            }
            .element-container:nth-of-type(3) button:hover {
                background-color: #c0c0c0 !important;
            }
            .element-container:nth-of-type(3) button:active {
                background-color: #a9a9a9 !important;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown('<div class="synonyms-button-container">', unsafe_allow_html=True)
    if "error_message" not in st.session_state:
        st.session_state["error_message"] = None
    if st.button(
        label="🔍  Retrieve Synonyms for\n\n\n Targets", key="retrieve_synonyms"
    ):
        if not targets_input.strip():
            st.session_state["error_message"] = "⚠️ Please enter at least one target!"
        else:
            target_list = [
                line.strip().split(",") for line in targets_input.strip().split("\n")
            ]
            synonyms_dict = {}
            error_found = False
            for entry in target_list:
                if len(entry) != 2:
                    st.session_state["error_message"] = (
                        f"❌ Invalid format for input: {' '.join(entry)}. Use format 'TARGET_NAME, TARGET_TYPE'."
                    )
                    error_found = True
                    continue
                target_name, target_type = entry[0].strip(), entry[1].strip().lower()
                if target_type not in [
                    "protein",
                    "gene",
                    "chemical",
                    "receptor",
                    "pathway",
                ]:
                    st.session_state["error_message"] = (
                        f"❌ Unknown target type: {target_type}. Choose from protein, gene, chemical, receptor, or pathway."
                    )
                    error_found = True
                    continue
                synonyms = retriever.get_target_synonyms(target_name, target_type)
                filtered_target_synonyms = [syn for syn in synonyms if syn.strip()]
                filtered_target_synonyms = filtered_target_synonyms[
                    :num_synonyms_per_target
                ]
                synonyms_dict[target_name] = [target_name] + filtered_target_synonyms
            if not error_found:
                st.session_state["error_message"] = None
            st.session_state["synonyms_dict"] = synonyms_dict
            st.rerun()
    if st.session_state["error_message"]:
        st.warning(st.session_state["error_message"])

# Additional keywords (unchanged)
additional_keywords_input = st.text_area(
    "🔗 Enter Other Keywords (One Keyword per Line) (Optional)"
)
additional_keywords_list = [
    keyword.strip()
    for keyword in additional_keywords_input.split("\n")
    if keyword.strip()
]
additional_condition = (
    f"AND ({' OR '.join([f'{kw}[Title/Abstract]' for kw in additional_keywords_list])})"
    if additional_keywords_list
    else ""
)

# Sidebar sliders and inputs (unchanged)
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

# Article type selection (unchanged)
article_types = [
    "Adaptive Clinical Trial",
    "Address",
    "Autobiography",
    "Bibliography",
    "Biography",
    "Books and Documents",
    "Case Reports",
    "Classical Article",
    "Clinical Conference",
    "Clinical Study",
    "Clinical Trial",
    "Clinical Trial Protocol",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Clinical Trial, Veterinary",
    "Collected Work",
    "Comment",
    "Comparative Study",
    "Congress",
    "Consensus Development Conference",
    "Consensus Development Conference, NIH",
    "Controlled Clinical Trial",
    "Corrected and Republished Article",
    "Dataset",
    "Dictionary",
    "Directory",
    "Duplicate Publication",
    "Editorial",
    "Electronic Supplementary Materials",
    "English Abstract",
    "Equivalence Trial",
    "Evaluation Study",
    "Expression of Concern",
    "Festschrift",
    "Government Publication",
    "Guideline",
    "Historical Article",
    "Interactive Tutorial",
    "Interview",
    "Introductory Journal Article",
    "Lecture",
    "Legal Case",
    "Legislation",
    "Letter",
    "Meta-Analysis",
    "Multicenter Study",
    "News",
    "Newspaper Article",
    "Observational Study",
    "Observational Study, Veterinary",
    "Overall",
    "Patient Education Handout",
    "Periodical Index",
    "Personal Narrative",
    "Portrait",
    "Practice Guideline",
    "Pragmatic Clinical Trial",
    "Preprint",
    "Published Erratum",
    "Randomized Controlled Trial",
    "Randomized Controlled Trial, Veterinary",
    "Research Support, American Recovery and Reinvestment Act",
    "Research Support, N.I.H., Extramural",
    "Research Support, N.I.H., Intramural",
    "Research Support, Non-U.S. Gov't",
    "Research Support, U.S. Gov't, Non-P.H.S.",
    "Research Support, U.S. Gov't, P.H.S.",
    "Research Support, U.S. Gov't",
    "Retracted Publication",
    "Retraction of Publication",
    "Review",
    "Scientific Integrity Review",
    "Systematic Review",
    "Technical Report",
    "Twin Study",
    "Validation Study",
    "Video-Audio Media",
    "Webcast",
]
selected_types = st.sidebar.multiselect("Filter by Article Type", article_types)
article_type_query = (
    " OR ".join([f'"{atype}"[Publication Type]' for atype in selected_types])
    if selected_types
    else ""
)

# Sidebar footer (unchanged)
st.sidebar.markdown(
    """
    <hr>
    <div style="text-align: center;">
        <p style="font-size: 12px; color: gray;">
            © 2024 Edelweiss Connect - Developed by Asmaa A. Abdelwahab
        </p>
    </div>
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


# Utility functions (unchanged)
def display_summary(compounds, targets):
    st.markdown("### Summary of Your Selections")
    st.markdown(f"**Email Address:** `{email if email else 'Not Provided'}`")
    st.markdown(f"**NCBI API Key:** `{api_key if api_key else 'Not Provided'}`")
    compounds_str = ", ".join(compounds)
    st.markdown(
        f"**Compounds List:** `{compounds_str if compounds_str else 'Not Provided'}`"
    )
    targets_str = ", ".join(targets) if targets else "Not Provided"
    st.markdown(f"**Interaction Targets List:** `{targets_str}`")
    keywords_str = (
        ", ".join(additional_keywords_list)
        if additional_keywords_list
        else "Not Provided"
    )
    st.markdown(f"**Additional Keywords:** `{keywords_str}`")
    st.markdown(
        f"**Number of Articles Per Compound-Target Pair:** `{n_articles_per_pair}`"
    )
    st.markdown(f"**Year Range:** `{start_year} to {end_year}`")
    st.markdown("---")


def generate_summary(task, compounds, targets):
    """
    Generate a plain-text summary of the user's selections for inclusion in the email.

    Parameters:
    task (dict): The task dictionary containing search parameters.
    compounds (list): List of compound names.
    targets (list): List of target names.

    Returns:
    str: A formatted summary string.
    """
    summary_lines = ["Summary of Your Selections:\n"]

    summary_lines.append(
        f"Email Address: {task['email'] if task['email'] else 'Not Provided'}"
    )
    summary_lines.append(
        f"NCBI API Key: {task['api_key'] if task['api_key'] else 'Not Provided'}"
    )

    compounds_str = ", ".join(compounds) if compounds else "Not Provided"
    summary_lines.append(f"Compounds List: {compounds_str}")

    targets_str = ", ".join(targets) if targets else "Not Provided"
    summary_lines.append(f"Interaction Targets List: {targets_str}")

    keywords_str = (
        ", ".join(task["additional_keywords_list"])
        if task["additional_keywords_list"]
        else "Not Provided"
    )
    summary_lines.append(f"Additional Keywords: {keywords_str}")

    summary_lines.append(
        f"Number of Articles Per Compound-Target Pair: {task['n_articles_per_pair']}"
    )
    summary_lines.append(f"Year Range: {task['start_year']} to {task['end_year']}")

    return "\n".join(summary_lines)


@st.fragment
def display_download_button(all_articles_df):
    csv = all_articles_df.to_csv(index=False).encode("utf-8")
    if st.download_button(
        label="📥 Download Combined Articles CSV",
        data=csv,
        file_name="combined_pubmed_articles.csv",
        mime="text/csv",
    ):
        st.success("✔️ Combined articles saved to 'combined_pubmed_articles.csv'")


def get_key_by_value(dictionary, value):
    """
    Safely get the key for a given value, handling lists comparison.
    """
    for key, val in dictionary.items():
        try:
            # Check if both are lists first
            if isinstance(val, list) and isinstance(value, list):
                if sorted(map(str.lower, val)) == sorted(map(str.lower, value)):
                    return key
            # If neither is a list, compare directly
            elif not isinstance(val, list) and not isinstance(value, list):
                if val == value:
                    return key
            # If types don't match (e.g., one is a list, the other isn't), skip
        except Exception as e:
            logging.error(f"Error matching key by value: {e}")
            continue
    return None


def is_valid_json5(text):
    try:
        json5.loads(text)
        return True
    except ValueError:
        return False


# Main search section (modified for queue system)
if st.button("🚀 Launch Search", help="Click to Start PubMed Search"):
    st.markdown("---")
    if not email or "@" not in email or "." not in email.split("@")[-1]:
        st.error("Please enter a valid email address.")
    elif not compounds_input:
        st.error("Please fill out the compound field.")
    else:
        # Process compounds
        compounds_list = [
            compound.strip()
            for compound in compounds_input.strip().split("\n")
            if compound.strip()
        ]
        if compounds_list:
            resolved_compounds = {
                compound: resolve_compound_name(compound) for compound in compounds_list
            }
            st.session_state["resolved_compounds"] = resolved_compounds
            st.session_state["compounds_text"] = format_compounds_json(
                resolved_compounds
            )

        if is_valid_json5(compounds_input):
            compounds_dict = json5.loads(compounds_input)
            compounds = list(compounds_dict.keys())
        else:
            resolved_compounds_list = st.session_state["resolved_compounds"]
            compounds_dict = {
                original: [resolved]
                for original, resolved in resolved_compounds.items()
            }
            compounds = compounds_list

        # Process targets
        if targets_input and targets_input.strip():
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
            targets_dict = {}
            targets = []

        # Display summary
        display_summary(compounds, targets if targets else ["No targets specified"])

        # Create task dictionary
        task = {
            "compounds_dict": compounds_dict,
            "targets_dict": targets_dict,
            "additional_keywords_list": additional_keywords_list,
            "n_articles_per_pair": n_articles_per_pair,
            "start_year": start_year,
            "end_year": end_year,
            "article_type_query": article_type_query,
            "email": email,
            "api_key": api_key,
        }

    # Check queue size before adding
    queue_size_before = task_queue.qsize()
    task_queue.put(task)

    if queue_size_before == 0:
        st.success("✅ Your search request has been submitted and is being processed.")
        st.info("📧 You will receive the results via email within an hour.")
    else:
        st.success("📝 Your request has been added to the queue.")
        st.info(
            f"⏳ You'll receive the results via email as soon as possible, and no later than 24 hours."
        )
