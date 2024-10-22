# 🧬 PubMed Retriever: Investigate Compound Interactions

**PubMed Retriever** is a Streamlit-based application designed to help researchers explore compound-gene interactions through comprehensive literature searches on PubMed. The tool facilitates the retrieval and analysis of scientific articles by querying combinations of compounds and genes, while also allowing the inclusion of specific keywords to refine the search.

## Features

- **Flexible Input Options:** Enter multiple compounds and genes of interest. Additionally, specify keywords to refine the search.
- **Automated Synonym Matching:** Automatically finds the most common synonyms for each compound to broaden the search.
- **Advanced Filtering:** Search articles within a specified year range.
- **Comprehensive Data Handling:** Filters duplicate results and handles missing or inconsistent data.
- **CSV Export:** Export search results to CSV for further analysis.

## How It Works

The app leverages several key technologies and data sources to retrieve relevant PubMed articles:

- **PubMed Integration:** Uses the NCBI PubMed API to fetch article metadata based on input criteria.
- **PubChem Synonym Finder:** Retrieves common synonyms for the input compounds from PubChem to expand the search.
- **Streamlit:** Provides a clean and interactive user interface for input, configuration, and output visualization.
- **Data Processing:** Combines and filters articles using Pandas to ensure accuracy and avoid duplicate entries.

## Installation

### Option 1: Local Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/your-repository/pubmed-retriever.git
   cd pubmed-retriever
   ```

2. **Set Up Virtual Environment:**

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install Dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

4. **Run the App:**

   ```bash
   streamlit run app.py
   ```

### Option 2: Running with Docker

1. **Build the Docker Image:**

   ```bash
   docker build -t pubmed-retriever .
   ```

2. **Run the Docker Container:**

   ```bash
   docker run -p 8501:8501 pubmed-retriever
   ```

3. **Access the App:**  
   Open a browser and go to `http://localhost:8501` to view and interact with the application.

## Usage

1. **Start the Application:** Open a terminal and run `streamlit run app.py` locally or launch the Docker container.
2. **Input Section:** 
   - **Compounds:** Enter compounds, each on a separate line.
   - **Genes:** (Optional) Enter genes, each on a separate line.
   - **Additional Keywords:** (Optional) Provide keywords to refine the search.
3. **Configuration:** Use the sidebar to configure various settings like the number of top synonyms per compound, date range, and number of articles to retrieve.
4. **Launch the Search:** Click on **"🚀 Launch Search"** to start the retrieval process.
5. **View and Export Results:** After processing, results will be displayed in a tabular format. You can export individual compound results or the combined results to CSV.

## Screenshots

1. **Application Input Page:**
   ![Application Input Page](images/1.png)

2. **Search Results Page:**
   ![Search Results Page](images/2.png)

## Dependencies

- **Python 3.7+**
- **Streamlit**: Web application framework for UI
- **Requests**: HTTP library to make API calls
- **Pandas**: Data analysis and manipulation
- **BioPython**: To interact with the NCBI databases
- **Metapub**: Fetch and process PubMed articles

## Docker Setup

The Dockerfile provided in this repository allows you to build and run the application in a containerized environment. Follow the steps below to use Docker:

1. **Build the Docker Image:**

   ```bash
   docker build -t pubmed-retriever .
   ```

2. **Run the Docker Container:**

   ```bash
   docker run -p 8501:8501 pubmed-retriever
   ```

3. **Access the App:**  
   Go to `http://localhost:8501` in your web browser.

### Dockerfile Overview

The Dockerfile uses a slim Python 3.9 base image and performs the following steps:

- Sets up a working directory.
- Copies the application code to the container.
- Installs the required Python packages.
- Exposes port `8501` for the Streamlit app.
- Runs the app on container startup using the Streamlit command.

## Contribution

Feel free to fork the repository and submit pull requests. If you encounter any issues or have suggestions, please create an issue in the repository.

## License

This project is licensed under the GPL-3.0 License. See the [LICENSE](./LICENSE) file for details.

## Acknowledgements

Special thanks to the creators of [Streamlit](https://streamlit.io/), [PubMed API](https://pubmed.ncbi.nlm.nih.gov/), and the [BioPython](https://biopython.org/) community for providing valuable resources.

## Contact

For any queries or assistance, please reach out through the repository's issue tracker or contact [@asmaa-a-abdelwahab](https://github.com/asmaa-a-abdelwahab).
