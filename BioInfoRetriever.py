import requests
import logging
import xmltodict
import time


class SynonymRetriever:
    def __init__(self):
        self.session = requests.Session()

    def _fetch_data(self, url, retries=3, delay=5):
        """Helper function to fetch data with retries."""
        for attempt in range(retries):
            try:
                response = self.session.get(url, timeout=10)
                response.raise_for_status()
                return (
                    response.json()
                    if response.headers.get("Content-Type", "").startswith(
                        "application/json"
                    )
                    else response.text
                )
            except requests.RequestException as e:
                logging.error("Error fetching data from: " + url)
                logging.error(
                    "Attempt "
                    + str(attempt + 1)
                    + " of "
                    + str(retries)
                    + ": "
                    + str(e)
                )
                if attempt < retries - 1:
                    time.sleep(delay * (2**attempt))
        return None

    ### **1️⃣ Proteins & Genes (UniProt, ChEMBL, NCBI Gene, HGNC)**
    def get_uniprot_synonyms(self, protein_name):
        """Retrieve synonyms for proteins/genes from UniProt."""
        url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_name}&fields=protein_name,gene_names"
        data = self._fetch_data(url)
        if not data or "results" not in data:
            return []

        synonyms = []
        for entry in data["results"]:
            if "protein_name" in entry:
                synonyms.append(entry["protein_name"])
            if "gene_names" in entry:
                synonyms.extend(
                    entry["gene_names"].split(" ")
                )  # Split space-separated gene names
        return list(set(synonyms))

    def get_ncbi_gene_synonyms(self, gene_symbol):
        """Retrieve synonyms for genes from NCBI Gene database."""
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_symbol}[Gene Name]&retmode=json"
        data = self._fetch_data(url)
        if not data or "esearchresult" not in data:
            return []

        gene_ids = data["esearchresult"].get("idlist", [])
        if not gene_ids:
            return []

        gene_id = gene_ids[0]
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json"
        summary_data = self._fetch_data(summary_url)
        if not summary_data:
            return []
        time.sleep(0.5)  # Add a small delay to avoid rate limiting
        return (
            summary_data.get("result", {})
            .get(gene_id, {})
            .get("otheraliases", "")
            .split(", ")
        )

    def get_hgnc_synonyms(self, gene_symbol):
        """Retrieve gene synonyms from HGNC."""
        url = f"https://rest.genenames.org/fetch/symbol/{gene_symbol}"
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)

        if response.status_code != 200:
            logging.error("HGNC API error: " + str(response.status_code))
            return []

        data = response.json()
        docs = data.get("response", {}).get("docs", [])

        if not docs:  # Ensure docs is not empty
            logging.warning("No synonyms found for: " + gene_symbol)
            return []

        return docs[0].get("alias_symbol", [])

    ### **2️⃣ Small Molecules & Drugs (PubChem, ChEMBL, DrugBank)**
    def get_pubchem_synonyms(self, chemical_name):
        """Retrieve synonyms for chemicals from PubChem."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/synonyms/JSON"
        data = self._fetch_data(url)
        if not data:
            return []

        return (
            data.get("InformationList", {})
            .get("Information", [{}])[0]
            .get("Synonym", [])
        )

    def get_chembl_synonyms(self, compound_name):
        """Retrieve synonyms for small molecules from ChEMBL."""
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__icontains={compound_name}"
        data = self._fetch_data(url)
        if not data or "molecules" not in data:
            return []

        return [
            compound["pref_name"]
            for compound in data["molecules"]
            if "pref_name" in compound
        ]

    ### **3️⃣ Transporters & Receptors (ChEMBL, UniProt)**
    def get_receptor_synonyms(self, receptor_name):
        """Retrieve synonyms for receptors from ChEMBL."""
        url = f"https://www.ebi.ac.uk/chembl/api/data/target.json?pref_name__icontains={receptor_name}"
        data = self._fetch_data(url)
        if not data or "targets" not in data:
            return []

        return [
            target["pref_name"] for target in data["targets"] if "pref_name" in target
        ]

    ### **4️⃣ Pathways (KEGG, Reactome, BioCyc)**
    def get_kegg_pathway_synonyms(self, pathway_name):
        """Retrieve synonyms for pathways from KEGG."""
        url = f"http://rest.kegg.jp/find/pathway/{pathway_name}"
        data = self._fetch_data(url)
        if not data:
            return []

        return [line.split("\t")[1] for line in data.split("\n") if "\t" in line]

    ### **Master Function to Get All Synonyms**
    def get_target_synonyms(self, entity_name, entity_type):
        """Retrieve synonyms from appropriate databases based on entity type."""
        synonyms = set()

        if entity_type in ["protein", "gene"]:
            synonyms.update(self.get_uniprot_synonyms(entity_name))
            synonyms.update(self.get_ncbi_gene_synonyms(entity_name))
            synonyms.update(self.get_hgnc_synonyms(entity_name))

        elif entity_type == "chemical":
            synonyms.update(self.get_pubchem_synonyms(entity_name))
            synonyms.update(self.get_chembl_synonyms(entity_name))

        elif entity_type == "receptor":
            synonyms.update(self.get_receptor_synonyms(entity_name))

        elif entity_type == "pathway":
            synonyms.update(self.get_kegg_pathway_synonyms(entity_name))

        return list(synonyms)


# # Example Usage
# # Initialize the synonym retriever
# retriever = SynonymRetriever()

# # 1️⃣ **Proteins & Genes**
# print("\n🔬 Testing Proteins & Genes")
# protein_synonyms = retriever.get_all_synonyms("BRAF", "gene")  # Example: BRAF gene
# print("BRAF Gene Synonyms:", protein_synonyms)

# protein_synonyms_2 = retriever.get_all_synonyms("TP53", "gene")  # Example: TP53 gene
# print("TP53 Gene Synonyms:", protein_synonyms_2)

# protein_synonyms_3 = retriever.get_all_synonyms(
#     "EGFR", "gene"
# )  # Epidermal Growth Factor Receptor
# print("EGFR Gene Synonyms:", protein_synonyms_3)

# protein_synonyms_4 = retriever.get_all_synonyms("MYC", "gene")  # MYC proto-oncogene
# print("MYC Gene Synonyms:", protein_synonyms_4)

# protein_synonyms_5 = retriever.get_all_synonyms(
#     "BRCA1", "gene"
# )  # BRCA1 (breast cancer gene)
# print("BRCA1 Gene Synonyms:", protein_synonyms_5)

# protein_synonyms_6 = retriever.get_all_synonyms(
#     "CDK2", "gene"
# )  # Cyclin-dependent kinase 2
# print("CDK2 Gene Synonyms:", protein_synonyms_6)

# protein_synonyms_7 = retriever.get_all_synonyms("KRAS", "gene")  # KRAS oncogene
# print("KRAS Gene Synonyms:", protein_synonyms_7)


# # 2️⃣ **Small Molecules & Drugs**
# print("\n🧪 Testing Small Molecules & Drugs")
# chemical_synonyms = retriever.get_all_synonyms(
#     "Aspirin", "chemical"
# )  # Example: Aspirin
# print("Aspirin Synonyms:", chemical_synonyms)

# chemical_synonyms_2 = retriever.get_all_synonyms(
#     "Ibuprofen", "chemical"
# )  # Example: Ibuprofen
# print("Ibuprofen Synonyms:", chemical_synonyms_2)

# chemical_synonyms_3 = retriever.get_all_synonyms(
#     "Paracetamol", "chemical"
# )  # Example: Paracetamol (Acetaminophen)
# print("Paracetamol Synonyms:", chemical_synonyms_3)

# chemical_synonyms_4 = retriever.get_all_synonyms(
#     "Metformin", "chemical"
# )  # Diabetes drug
# print("Metformin Synonyms:", chemical_synonyms_4)

# chemical_synonyms_5 = retriever.get_all_synonyms("Penicillin", "chemical")  # Antibiotic
# print("Penicillin Synonyms:", chemical_synonyms_5)

# chemical_synonyms_6 = retriever.get_all_synonyms("Caffeine", "chemical")  # Stimulant
# print("Caffeine Synonyms:", chemical_synonyms_6)

# chemical_synonyms_7 = retriever.get_all_synonyms(
#     "Atorvastatin", "chemical"
# )  # Cholesterol-lowering drug
# print("Atorvastatin Synonyms:", chemical_synonyms_7)


# # 3️⃣ **Transporters & Receptors**
# print("\n⚡ Testing Transporters & Receptors")
# receptor_synonyms = retriever.get_all_synonyms(
#     "Dopamine receptor", "receptor"
# )  # Example: Dopamine receptor
# print("Dopamine Receptor Synonyms:", receptor_synonyms)

# receptor_synonyms_2 = retriever.get_all_synonyms(
#     "Serotonin transporter", "receptor"
# )  # Example: Serotonin transporter
# print("Serotonin Transporter Synonyms:", receptor_synonyms_2)

# receptor_synonyms_3 = retriever.get_all_synonyms(
#     "Histamine receptor", "receptor"
# )  # Histamine receptor
# print("Histamine Receptor Synonyms:", receptor_synonyms_3)

# receptor_synonyms_4 = retriever.get_all_synonyms(
#     "Insulin receptor", "receptor"
# )  # Insulin receptor
# print("Insulin Receptor Synonyms:", receptor_synonyms_4)

# receptor_synonyms_5 = retriever.get_all_synonyms(
#     "Glutamate receptor", "receptor"
# )  # Glutamate receptor
# print("Glutamate Receptor Synonyms:", receptor_synonyms_5)

# receptor_synonyms_6 = retriever.get_all_synonyms(
#     "GABA receptor", "receptor"
# )  # GABA receptor
# print("GABA Receptor Synonyms:", receptor_synonyms_6)

# receptor_synonyms_7 = retriever.get_all_synonyms(
#     "Nicotinic acetylcholine receptor", "receptor"
# )  # Nicotinic acetylcholine receptor
# print("Nicotinic Acetylcholine Receptor Synonyms:", receptor_synonyms_7)


# # 4️⃣ **Pathways**
# print("\n🧬 Testing Pathways")
# pathway_synonyms = retriever.get_all_synonyms(
#     "Glycolysis", "pathway"
# )  # Example: Glycolysis pathway
# print("Glycolysis Pathway Synonyms:", pathway_synonyms)

# pathway_synonyms_2 = retriever.get_all_synonyms(
#     "Apoptosis", "pathway"
# )  # Example: Apoptosis pathway
# print("Apoptosis Pathway Synonyms:", pathway_synonyms_2)

# pathway_synonyms_3 = retriever.get_all_synonyms(
#     "TCA Cycle", "pathway"
# )  # Tricarboxylic Acid Cycle (Krebs Cycle)
# print("TCA Cycle Pathway Synonyms:", pathway_synonyms_3)

# pathway_synonyms_4 = retriever.get_all_synonyms(
#     "MAPK signaling pathway", "pathway"
# )  # Mitogen-activated protein kinase pathway
# print("MAPK Signaling Pathway Synonyms:", pathway_synonyms_4)

# pathway_synonyms_5 = retriever.get_all_synonyms(
#     "mTOR signaling pathway", "pathway"
# )  # mTOR signaling
# print("mTOR Signaling Pathway Synonyms:", pathway_synonyms_5)

# pathway_synonyms_6 = retriever.get_all_synonyms(
#     "Wnt signaling pathway", "pathway"
# )  # Wnt signaling
# print("Wnt Signaling Pathway Synonyms:", pathway_synonyms_6)

# pathway_synonyms_7 = retriever.get_all_synonyms(
#     "Oxidative phosphorylation", "pathway"
# )  # Energy production pathway
# print("Oxidative Phosphorylation Pathway Synonyms:", pathway_synonyms_7)
