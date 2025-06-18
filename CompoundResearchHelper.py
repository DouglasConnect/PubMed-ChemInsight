import logging
import time
import requests
import pandas as pd
from metapub import PubMedFetcher


class CompoundResearchHelper:
    """
    A professional class to fetch compound synonyms and retrieve PubMed articles.
    Ensures strict input validation, optimized query handling, and structured logging.
    """

    def __init__(self, retmax: int = 1000, api_key: str = None) -> None:
        """
        Initialize the helper with PubMedFetcher and settings.

        Args:
            retmax (int): Maximum number of articles per query.
            api_key (str): Optional NCBI API Key for increased rate limits.
        """
        self.pubmed = PubMedFetcher(api_key=api_key) if api_key else PubMedFetcher()
        self.retmax = retmax
        self.articleList = []
        self.pubchem_base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        logging.info(f"CompoundResearchHelper initialized with retmax={retmax}")

    def _fetch_data(self, url: str) -> dict:
        """Fetch JSON data from a URL with error handling."""
        if not isinstance(url, str):
            logging.error("URL must be a string.")
            return {}
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as err:
            logging.error(f"Request error: {url}: {err}")
            return {}

    def _clean_text(self, text: str) -> str:
        """Remove unwanted characters from search terms."""
        return "".join(c for c in text if c.isalnum() or c.isspace())

    def get_pubchem_synonyms(self, chemical_name: str) -> list:
        """Retrieve synonyms from PubChem."""
        if not isinstance(chemical_name, str):
            logging.error("Chemical name must be a string.")
            return []
        try:
            cid_url = f"{self.pubchem_base_url}/compound/name/{chemical_name}/cids/JSON"
            cid_data = self._fetch_data(cid_url)
            cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]
            if not cid:
                logging.warning(f"No CID found for {chemical_name}.")
                return []
            synonyms_url = f"{self.pubchem_base_url}/compound/cid/{cid}/synonyms/JSON"
            synonyms_data = self._fetch_data(synonyms_url)
            synonyms = (
                synonyms_data.get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )
            logging.info(f"Found {len(synonyms)} synonyms for {chemical_name}.")
            return synonyms
        except Exception as e:
            logging.error(f"PubChem synonym retrieval error for {chemical_name}: {e}")
            return []

    def get_compound_synonyms(self, compound_name: str) -> list:
        """Retrieve synonyms plus cleaned original name."""
        synonyms = self.get_pubchem_synonyms(compound_name)
        all_names = [self._clean_text(compound_name)] + [
            self._clean_text(s) for s in synonyms
        ]
        unique_synonyms = list(set(name for name in all_names if name))
        return unique_synonyms if unique_synonyms else [compound_name]

    def fetch_articles(
        self,
        search_term: str,
        retmax: int = 1000,
        start_year: int = 2000,
        end_year: int = None,
        email: str = None,
        api_key: str = None,
        retries: int = 3,
        backoff_factor: float = 2.0,
    ) -> list:
        if not isinstance(search_term, str):
            logging.error("Search term must be a string.")
            return []

        date_range = (
            f' AND ("{start_year}/01/01"[PDat] : "{end_year}/12/31"[PDat])'
            if end_year
            else ""
        )
        full_search_term = f"{search_term}{date_range}"
        logging.debug(f"Full search term: {full_search_term}")

        # Retry fetching PMIDs
        for attempt in range(retries):
            try:
                pmids = self.pubmed.pmids_for_query(
                    full_search_term, retmax=retmax, sort="relevance"
                )
                if not pmids:
                    logging.warning(
                        f"No PMIDs found for search term: {full_search_term}"
                    )
                else:
                    logging.info(
                        f"Fetched {len(pmids)} PMIDs for search: {search_term}"
                    )
                break
            except Exception as e:
                if attempt < retries - 1:
                    wait_time = backoff_factor**attempt
                    logging.warning(
                        f"Error fetching PMIDs for '{search_term}': {e}. Retrying in {wait_time} seconds..."
                    )
                    time.sleep(wait_time)
                else:
                    logging.error(
                        f"Failed to fetch PMIDs after {retries} attempts: {e}"
                    )
                    return []

        articles = []
        for pmid in pmids:
            try:
                article = self.pubmed.article_by_pmid(pmid)
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
                article_dict = {k: getattr(article, k, None) for k in keys_to_keep}
                # Normalize the year field to an integer or None
                year_value = article_dict["year"]
                if year_value is None or year_value == "":
                    article_dict["year"] = None
                elif isinstance(year_value, int):
                    article_dict["year"] = year_value
                else:
                    import re

                    year_str = str(year_value)
                    year_match = re.match(r"^\d{4}", year_str)
                    if year_match:
                        try:
                            article_dict["year"] = int(year_match.group(0))
                        except (ValueError, TypeError):
                            logging.warning(
                                f"Cannot convert year for PMID {pmid}: {year_value!r}"
                            )
                            article_dict["year"] = None
                    else:
                        logging.warning(
                            f"Invalid year format for PMID {pmid}: {year_value!r}"
                        )
                        article_dict["year"] = None
                articles.append(article_dict)
            except Exception as e:
                logging.warning(f"Failed to process PMID {pmid}: {e}")

        if not articles and pmids:
            logging.warning(
                f"Fetched {len(pmids)} PMIDs but no articles were processed successfully."
            )
        return articles

    def select_top_articles(self, df: pd.DataFrame, n_articles: int) -> pd.DataFrame:
        """Select top recent articles."""
        if not isinstance(df, pd.DataFrame):
            logging.error("Input must be a pandas DataFrame.")
            return pd.DataFrame()
        if "year" in df.columns:
            df["year"] = pd.to_numeric(df["year"], errors="coerce")
        df = df.drop_duplicates(subset=["pmid"]).sort_values(by="year", ascending=False)
        top_articles = df.head(n_articles) if n_articles > 0 else pd.DataFrame()
        logging.info(f"Selected {len(top_articles)} top articles.")
        return top_articles

    ######################################################################################
    def _escape_pubmed_query(self, term: str) -> str:
        """Escape special characters for PubMed query syntax."""
        import re

        return re.sub(r"([^\w\s])", r"\\\1", term)

    def build_term_query(self, term: str, fields: list[str]) -> str:
        """Build OR-separated field query for a single term."""
        clean = self._escape_pubmed_query(self._clean_text(term))
        return " OR ".join(f'"{clean}"[{field}]' for field in fields)

    def get_batch_query(self, terms: list[str], fields: list[str]) -> str:
        """
        Construct a PubMed query for a batch of terms without synonyms.

        Args:
            terms (list[str]): Terms to include in the query.
            fields (list[str]): PubMed fields like 'Title/Abstract', 'MeSH Terms', etc.

        Returns:
            str: Combined query string for PubMed.
        """
        queries = []
        for term in terms:
            clean = self._escape_pubmed_query(self._clean_text(term))
            queries.extend(f'"{clean}"[{field}]' for field in fields)
        return " OR ".join(queries)

    def _format_article_type_filter(self, types: list[str]) -> str:
        allowed = {"review", "clinical trial", "case reports"}
        valid = [t for t in types if t.lower() in allowed]
        return (
            f" AND ({' OR '.join(f'{t}[Publication Type]' for t in valid)})"
            if valid
            else ""
        )

    def process_compound_and_targets(
        self,
        compounds: list,
        genes: list,
        start_year: int,
        end_year: int,
        additional_condition: str,
        n_articles: int,
        article_type_query: str = None,
        fields: list = ["Title/Abstract", "MeSH Terms", "Substance Name"],
    ) -> pd.DataFrame:
        """Main function to search PubMed and retrieve articles."""
        if not (
            isinstance(compounds, list) and all(isinstance(c, str) for c in compounds)
        ):
            logging.error("Compounds must be a list of strings.")
            return pd.DataFrame()
        if genes and not (
            isinstance(genes, list) and all(isinstance(g, str) for g in genes)
        ):
            logging.error("Genes must be a list of strings.")
            return pd.DataFrame()
        if start_year > end_year:
            logging.error("Start year cannot be after end year.")
            return pd.DataFrame()

        self.articleList = []
        batch_size = 5
        article_type_condition = (
            self._format_article_type_filter([article_type_query])
            if article_type_query
            else ""
        )

        queries = []
        for c_batch in (
            compounds[i : i + batch_size] for i in range(0, len(compounds), batch_size)
        ):
            compound_query = self.get_batch_query(c_batch, fields)
            if genes:
                for g_batch in (
                    genes[i : i + batch_size] for i in range(0, len(genes), batch_size)
                ):
                    gene_query = self.get_batch_query(g_batch, ["Title/Abstract"])
                    full_query = f"(({compound_query}) AND ({gene_query})){additional_condition}{article_type_condition}"
                    queries.append(full_query)
            else:
                full_query = (
                    f"({compound_query}){additional_condition}{article_type_condition}"
                )
                queries.append(full_query)

        for idx, query in enumerate(queries, 1):
            self.articleList.extend(
                self.fetch_articles(
                    query, retmax=self.retmax, start_year=start_year, end_year=end_year
                )
            )
            if idx % 3 == 0:
                time.sleep(1)

        if self.articleList:
            df = pd.DataFrame(self.articleList)
            return self.select_top_articles(df, n_articles)

        logging.warning("No articles retrieved.")
        return pd.DataFrame()

    ######################################################################################

    # def process_compound_and_targets(
    #     self,
    #     compounds: list,
    #     genes: list,
    #     start_year: int,
    #     end_year: int,
    #     additional_condition: str,
    #     n_articles: int,
    #     article_type_query: str = None,
    # ) -> pd.DataFrame:
    #     """Main function to search PubMed and retrieve articles."""
    #     if not (
    #         isinstance(compounds, list) and all(isinstance(c, str) for c in compounds)
    #     ):
    #         logging.error("Compounds must be a list of strings.")
    #         return pd.DataFrame()
    #     if genes and not (
    #         isinstance(genes, list) and all(isinstance(g, str) for g in genes)
    #     ):
    #         logging.error("Genes must be a list of strings.")
    #         return pd.DataFrame()
    #     if start_year > end_year:
    #         logging.error("Start year cannot be after end year.")
    #         return pd.DataFrame()

    #     self.articleList = []
    #     batch_size = 5
    #     queries = []
    #     article_type_condition = (
    #         f" AND ({article_type_query})" if article_type_query else ""
    #     )

    #     if genes:
    #         for c_batch in (
    #             compounds[i : i + batch_size]
    #             for i in range(0, len(compounds), batch_size)
    #         ):
    #             for g_batch in (
    #                 genes[i : i + batch_size] for i in range(0, len(genes), batch_size)
    #             ):
    #                 query = f"(({' OR '.join(f'{self._clean_text(c)}[Title/Abstract]' for c in c_batch)}) AND ({' OR '.join(f'{self._clean_text(g)}[Title/Abstract]' for g in g_batch)})){additional_condition}{article_type_condition}"
    #                 queries.append(query)
    #     else:
    #         for c_batch in (
    #             compounds[i : i + batch_size]
    #             for i in range(0, len(compounds), batch_size)
    #         ):
    #             query = f"({' OR '.join(f'{self._clean_text(c)}[Title/Abstract]' for c in c_batch)}){additional_condition}{article_type_condition}"
    #             queries.append(query)

    #     for idx, query in enumerate(queries, 1):
    #         self.articleList.extend(
    #             self.fetch_articles(
    #                 query, retmax=self.retmax, start_year=start_year, end_year=end_year
    #             )
    #         )
    #         if idx % 3 == 0:
    #             time.sleep(1)

    #     if self.articleList:
    #         df = pd.DataFrame(self.articleList)
    #         return self.select_top_articles(df, n_articles)

    #     logging.warning("No articles retrieved.")
    #     return pd.DataFrame()


# import logging
# import time
# import requests
# import pandas as pd
# from itertools import product
# from metapub import PubMedFetcher


# class CompoundResearchHelper:
#     """A class to fetch compound synonyms and retrieve relevant articles from PubMed."""

#     def __init__(self, retmax=1000):
#         """
#         Initialize the CompoundResearchHelper instance with default settings.

#         Parameters:
#         retmax (int): The maximum number of articles to retrieve per query. Defaults to 1000.

#         Attributes:
#         pubmed (PubMedFetcher): An instance of PubMedFetcher for retrieving PubMed articles.
#         retmax (int): Stores the maximum number of articles to retrieve.
#         articleList (list): A list to store fetched articles.
#         """

#         self.pubmed = PubMedFetcher()
#         self.retmax = retmax
#         self.articleList = []

#         self.pubchem_base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
#         self.chembl_base_url = "https://www.ebi.ac.uk/chembl/api/data"

#     def _fetch_data(self, url):
#         """Helper function to handle API requests."""
#         try:
#             response = requests.get(url)
#             response.raise_for_status()
#             return response.json()
#         except requests.RequestException as err:
#             logging.error(f"Error fetching data from {url}: {err}")
#             return {}

#     def get_pubchem_synonyms(self, chemical_name):
#         """Retrieve synonyms for chemicals from PubChem."""
#         try:
#             # Step 1: Get CID
#             cid_url = f"{self.pubchem_base_url}/compound/name/{chemical_name}/cids/JSON"
#             cid_data = self._fetch_data(cid_url)
#             cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]
#             if not cid:
#                 logging.warning(f"No CID found for {chemical_name} in PubChem.")
#                 return []

#             # Step 2: Retrieve synonyms using CID
#             synonyms_url = f"{self.pubchem_base_url}/compound/cid/{cid}/synonyms/JSON"
#             synonyms_data = self._fetch_data(synonyms_url)

#             return (
#                 synonyms_data.get("InformationList", {})
#                 .get("Information", [{}])[0]
#                 .get("Synonym", [])
#             )
#         except Exception as e:
#             logging.error(f"Error retrieving PubChem synonyms for {chemical_name}: {e}")
#             return []

#     def get_chembl_synonyms(self, compound_name):
#         """Retrieve synonyms for small molecules from ChEMBL."""
#         url = (
#             f"{self.chembl_base_url}/molecule.json?pref_name__icontains={compound_name}"
#         )
#         data = self._fetch_data(url)
#         if not data or "molecules" not in data:
#             logging.warning(f"No ChEMBL data found for {compound_name}.")
#             return []

#         synonyms = [
#             compound["pref_name"]
#             for compound in data["molecules"]
#             if "pref_name" in compound
#         ]
#         return list(set(synonyms))  # Ensure unique synonyms

#     def get_compound_synonyms(self, compound_name):
#         """
#         Retrieves synonyms for the given compound name from PubChem and ChEMBL.

#         Parameters:
#         compound_name (str): The name of the compound to retrieve synonyms for.

#         Returns:
#         list: A combined list of unique synonyms from PubChem and ChEMBL.
#         """
#         pubchem_synonyms = self.get_pubchem_synonyms(compound_name)
#         # chembl_synonyms = self.get_chembl_synonyms(compound_name)

#         # Combine and remove duplicates
#         all_synonyms = list(set([compound_name] + pubchem_synonyms))
#         if all_synonyms:
#             return all_synonyms
#         else:
#             [compound_name]

#     def fetch_articles(self, search_term, retmax=1000, start_year=2000, end_year=None):
#         """
#         Fetches PubMed articles for a given search term, with optional filters for date range.

#         Parameters:
#         search_term (str): The search term to query PubMed with.
#         retmax (int): The maximum number of results to return from the query (default: 1000).
#         start_year (int): The start of the year range to filter the results by (default: 2000).
#         end_year (int): The end of the year range to filter the results by (default: None, i.e., current year).

#         Returns:
#         list: A list of dictionaries, each containing information about a single article.
#         """

#         date_range = (
#             f' AND ("{start_year}/01/01"[PDat] : "{end_year}/12/31"[PDat])'
#             if end_year
#             else ""
#         )
#         full_search_term = f"{search_term}{date_range}"

#         try:
#             pmids = self.pubmed.pmids_for_query(
#                 full_search_term, retmax=retmax, sort="relevance"
#             )
#         except Exception as e:
#             logging.error(f"Error fetching pmids for {search_term}: {e}")
#             return []

#         articles = []
#         for pmid in pmids:
#             try:
#                 article = self.pubmed.article_by_pmid(pmid)
#                 article_dict = article.to_dict()
#                 keys_to_keep = [
#                     "title",
#                     "pmid",
#                     "url",
#                     "authors",
#                     "doi",
#                     "pmc",
#                     "issn",
#                     "mesh",
#                     "chemicals",
#                     "journal",
#                     "abstract",
#                     "year",
#                     "publication_types",
#                 ]
#                 article_dict = {k: article_dict.get(k, None) for k in keys_to_keep}
#                 articles.append(article_dict)
#             except Exception as e:
#                 logging.error(f"Error processing article with pmid {pmid}: {e}")
#         return articles

#     def filter_articles_by_recent(self, df, top_n, bottom_n):
#         """
#         Filters and returns a DataFrame containing the top and bottom recent articles based on the year.

#         Parameters:
#         df (pd.DataFrame): DataFrame containing article data with a 'year' column.
#         top_n (int): Number of top recent articles to select based on the year.
#         bottom_n (int): Number of bottom recent articles to select based on the year.

#         Returns:
#         pd.DataFrame: A DataFrame containing the top and bottom recent articles, deduplicated by 'title' and 'pmid'.
#         """

#         if "year" in df.columns:
#             df["year"] = pd.to_numeric(df["year"], errors="coerce")

#         df_sorted = df.sort_values(by="year", ascending=False)

#         df_sorted = df_sorted.map(
#             lambda x: str(x) if isinstance(x, (list, dict)) else x
#         )

#         top_recent = df_sorted.head(top_n) if top_n > 0 else pd.DataFrame()
#         bottom_recent = df_sorted.tail(bottom_n) if bottom_n > 0 else pd.DataFrame()

#         filtered_df = pd.concat([top_recent, bottom_recent]).drop_duplicates(
#             subset=["title", "pmid"]
#         )

#         return filtered_df

#     def process_compound_and_targets(
#         self,
#         compounds,
#         genes,
#         start_year,
#         end_year,
#         additional_condition,
#         top_n,
#         bottom_n,
#         article_type_query=None,
#     ):
#         """
#         Process a compound and optional list of genes by performing a PubMed search
#         with the compound and genes, and then filtering the results to include only
#         the top and bottom recent articles.

#         Parameters:
#         compound (str): The name of the compound to search for.
#         genes (list): A list of gene symbols to search for in conjunction with the compound.
#         start_year (int): The starting year to use in the search.
#         end_year (int): The ending year to use in the search.
#         additional_condition (str): An additional condition to use in the search.
#         top_n (int): The number of top recent articles to select based on the year.
#         bottom_n (int): The number of bottom recent articles to select based on the year.
#         article_type_query (str, optional): PubMed article type filter in query format.

#         Returns:
#         pd.DataFrame: A DataFrame containing the top and bottom recent articles, deduplicated by 'title' and 'pmid'.
#         """
#         logging.info(f"Processing compound: {compounds[0]}")

#         batch_size = 5  # Adjust batch size to avoid long queries
#         queries = []

#         # Append article type query if provided
#         article_type_condition = (
#             f" AND ({article_type_query})" if article_type_query else ""
#         )

#         if genes:
#             for compound_chunk in [
#                 compounds[i : i + batch_size]
#                 for i in range(0, len(compounds), batch_size)
#             ]:
#                 for gene_chunk in [
#                     genes[i : i + batch_size] for i in range(0, len(genes), batch_size)
#                 ]:
#                     query = (
#                         f"(({' OR '.join([f'{compound}[Title/Abstract]' for compound in compound_chunk])}) AND "
#                         f"({' OR '.join([f'{gene}[Title/Abstract]' for gene in gene_chunk])})) "
#                         f"{additional_condition}{article_type_condition}"
#                     )
#                     queries.append(query)
#         else:
#             for compound_chunk in [
#                 compounds[i : i + batch_size]
#                 for i in range(0, len(compounds), batch_size)
#             ]:
#                 query = (
#                     f"({' OR '.join([f'{compound}[Title/Abstract]' for compound in compound_chunk])}) "
#                     f"{additional_condition}{article_type_condition}"
#                 )
#                 queries.append(query)

#         query_count = 0
#         for query in queries:
#             self.articleList.extend(
#                 self.fetch_articles(
#                     query, retmax=self.retmax, start_year=start_year, end_year=end_year
#                 )
#             )

#             # Increment the query counter
#             query_count += 1

#             # Check if 3 queries have been made
#             if query_count == 3:
#                 # Sleep for the desired amount of time
#                 time.sleep(1)  # Adjust the sleep time (e.g., 1 second) as needed

#                 # Reset the counter
#                 query_count = 0

#         if self.articleList:
#             df = pd.DataFrame(self.articleList)
#             df_filtered = self.filter_articles_by_recent(df, top_n, bottom_n)
#             return df_filtered
#         return pd.DataFrame()
