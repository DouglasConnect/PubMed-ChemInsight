import logging
import time
import requests
import pandas as pd
from itertools import product
from metapub import PubMedFetcher


class CompoundResearchHelper:
    """A class to fetch compound synonyms and retrieve relevant articles from PubMed."""

    def __init__(self, retmax=1000):
        """
        Initialize the CompoundResearchHelper instance with default settings.

        Parameters:
        retmax (int): The maximum number of articles to retrieve per query. Defaults to 1000.

        Attributes:
        pubmed (PubMedFetcher): An instance of PubMedFetcher for retrieving PubMed articles.
        retmax (int): Stores the maximum number of articles to retrieve.
        articleList (list): A list to store fetched articles.
        """

        self.pubmed = PubMedFetcher()
        self.retmax = retmax
        self.articleList = []

        self.pubchem_base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.chembl_base_url = "https://www.ebi.ac.uk/chembl/api/data"

    def _fetch_data(self, url):
        """Helper function to handle API requests."""
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as err:
            logging.error(f"Error fetching data from {url}: {err}")
            return {}

    def get_pubchem_synonyms(self, chemical_name):
        """Retrieve synonyms for chemicals from PubChem."""
        try:
            # Step 1: Get CID
            cid_url = f"{self.pubchem_base_url}/compound/name/{chemical_name}/cids/JSON"
            cid_data = self._fetch_data(cid_url)
            cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]
            if not cid:
                logging.warning(f"No CID found for {chemical_name} in PubChem.")
                return []

            # Step 2: Retrieve synonyms using CID
            synonyms_url = f"{self.pubchem_base_url}/compound/cid/{cid}/synonyms/JSON"
            synonyms_data = self._fetch_data(synonyms_url)

            return (
                synonyms_data.get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )
        except Exception as e:
            logging.error(f"Error retrieving PubChem synonyms for {chemical_name}: {e}")
            return []

    def get_chembl_synonyms(self, compound_name):
        """Retrieve synonyms for small molecules from ChEMBL."""
        url = (
            f"{self.chembl_base_url}/molecule.json?pref_name__icontains={compound_name}"
        )
        data = self._fetch_data(url)
        if not data or "molecules" not in data:
            logging.warning(f"No ChEMBL data found for {compound_name}.")
            return []

        synonyms = [
            compound["pref_name"]
            for compound in data["molecules"]
            if "pref_name" in compound
        ]
        return list(set(synonyms))  # Ensure unique synonyms

    def get_compound_synonyms(self, compound_name):
        """
        Retrieves synonyms for the given compound name from PubChem and ChEMBL.

        Parameters:
        compound_name (str): The name of the compound to retrieve synonyms for.

        Returns:
        list: A combined list of unique synonyms from PubChem and ChEMBL.
        """
        pubchem_synonyms = self.get_pubchem_synonyms(compound_name)
        # chembl_synonyms = self.get_chembl_synonyms(compound_name)

        # Combine and remove duplicates
        all_synonyms = list(set([compound_name] + pubchem_synonyms))
        if all_synonyms:
            return all_synonyms
        else:
            [compound_name]

    def fetch_articles(self, search_term, retmax=1000, start_year=2000, end_year=None):
        """
        Fetches PubMed articles for a given search term, with optional filters for date range.

        Parameters:
        search_term (str): The search term to query PubMed with.
        retmax (int): The maximum number of results to return from the query (default: 1000).
        start_year (int): The start of the year range to filter the results by (default: 2000).
        end_year (int): The end of the year range to filter the results by (default: None, i.e., current year).

        Returns:
        list: A list of dictionaries, each containing information about a single article.
        """

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
                    "publication_types",
                ]
                article_dict = {k: article_dict.get(k, None) for k in keys_to_keep}
                articles.append(article_dict)
            except Exception as e:
                logging.error(f"Error processing article with pmid {pmid}: {e}")
        return articles

    def filter_articles_by_recent(self, df, top_n, bottom_n):
        """
        Filters and returns a DataFrame containing the top and bottom recent articles based on the year.

        Parameters:
        df (pd.DataFrame): DataFrame containing article data with a 'year' column.
        top_n (int): Number of top recent articles to select based on the year.
        bottom_n (int): Number of bottom recent articles to select based on the year.

        Returns:
        pd.DataFrame: A DataFrame containing the top and bottom recent articles, deduplicated by 'title' and 'pmid'.
        """

        if "year" in df.columns:
            df["year"] = pd.to_numeric(df["year"], errors="coerce")

        df_sorted = df.sort_values(by="year", ascending=False)

        df_sorted = df_sorted.map(
            lambda x: str(x) if isinstance(x, (list, dict)) else x
        )

        top_recent = df_sorted.head(top_n) if top_n > 0 else pd.DataFrame()
        bottom_recent = df_sorted.tail(bottom_n) if bottom_n > 0 else pd.DataFrame()

        filtered_df = pd.concat([top_recent, bottom_recent]).drop_duplicates(
            subset=["title", "pmid"]
        )

        return filtered_df

    def process_compound_and_targets(
        self,
        compounds,
        genes,
        start_year,
        end_year,
        additional_condition,
        top_n,
        bottom_n,
        article_type_query=None,
    ):
        """
        Process a compound and optional list of genes by performing a PubMed search
        with the compound and genes, and then filtering the results to include only
        the top and bottom recent articles.

        Parameters:
        compound (str): The name of the compound to search for.
        genes (list): A list of gene symbols to search for in conjunction with the compound.
        start_year (int): The starting year to use in the search.
        end_year (int): The ending year to use in the search.
        additional_condition (str): An additional condition to use in the search.
        top_n (int): The number of top recent articles to select based on the year.
        bottom_n (int): The number of bottom recent articles to select based on the year.
        article_type_query (str, optional): PubMed article type filter in query format.

        Returns:
        pd.DataFrame: A DataFrame containing the top and bottom recent articles, deduplicated by 'title' and 'pmid'.
        """
        logging.info(f"Processing compound: {compounds[0]}")

        batch_size = 5  # Adjust batch size to avoid long queries
        queries = []

        # Append article type query if provided
        article_type_condition = (
            f" AND ({article_type_query})" if article_type_query else ""
        )

        if genes:
            for compound_chunk in [
                compounds[i : i + batch_size]
                for i in range(0, len(compounds), batch_size)
            ]:
                for gene_chunk in [
                    genes[i : i + batch_size] for i in range(0, len(genes), batch_size)
                ]:
                    query = (
                        f"(({' OR '.join([f'{compound}[Title/Abstract]' for compound in compound_chunk])}) AND "
                        f"({' OR '.join([f'{gene}[Title/Abstract]' for gene in gene_chunk])})) "
                        f"{additional_condition}{article_type_condition}"
                    )
                    queries.append(query)
        else:
            for compound_chunk in [
                compounds[i : i + batch_size]
                for i in range(0, len(compounds), batch_size)
            ]:
                query = (
                    f"({' OR '.join([f'{compound}[Title/Abstract]' for compound in compound_chunk])}) "
                    f"{additional_condition}{article_type_condition}"
                )
                queries.append(query)

        query_count = 0
        for query in queries:
            self.articleList.extend(
                self.fetch_articles(
                    query, retmax=self.retmax, start_year=start_year, end_year=end_year
                )
            )

            # Increment the query counter
            query_count += 1

            # Check if 3 queries have been made
            if query_count == 3:
                # Sleep for the desired amount of time
                time.sleep(1)  # Adjust the sleep time (e.g., 1 second) as needed

                # Reset the counter
                query_count = 0

        if self.articleList:
            df = pd.DataFrame(self.articleList)
            df_filtered = self.filter_articles_by_recent(df, top_n, bottom_n)
            return df_filtered
        return pd.DataFrame()
