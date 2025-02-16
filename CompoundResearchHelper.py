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

    def get_compound_synonyms(self, compound_name):
        """
        Retrieves synonyms for the given compound name from PubChem.

        Parameters:
        compound_name (str): The name of the compound to retrieve synonyms for.

        Returns:
        list: A list of synonyms for the given compound name. If no synonyms can be found, or if an error occurs, an empty list is returned.

        Raises:
        requests.exceptions.HTTPError: Raised if an HTTP error occurs during the request.
        Exception: Raised if any other error occurs during the request.
        """
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
        compound,
        genes,
        start_year,
        end_year,
        additional_condition,
        top_n,
        bottom_n,
        top_syn,
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
        top_syn (int): The number of top synonyms to use in the search.

        Returns:
        pd.DataFrame: A DataFrame containing the top and bottom recent articles, deduplicated by 'title' and 'pmid'.
        """
        logging.info(f"Processing compound: {compound}")

        top_synonyms = set([compound] + self.get_compound_synonyms(compound)[:top_syn])

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
