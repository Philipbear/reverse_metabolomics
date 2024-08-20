
This is to generate the library for reverse metabolomics in batch.

- Input:
  - **mzML** files
  - **csv** files (compound list)

- Output:
  - one library **tsv** file
  - `metadata` folder
    - entire metadata table for each mzML
    - EIC plots for each compound in each mzML
    - MS2 spectra (selected & unselected) for each mzML

### Run the workflow
1. Clone the repository to local.

2. Install the dependencies.
```bash
pip install -r requirements.txt
```
3. prepare all the files in `input` folder.

4. `run.py` is the main script to run the workflow.

