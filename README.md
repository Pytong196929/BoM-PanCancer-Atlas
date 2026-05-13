**A Pan-Cancer Single-Cell Atlas of Bone Metastasis**
An interactive web portal for exploring cell types, gene expression, and compositional landscapes across multiple cancer types with bone metastasis at single-cell resolution.

BoM-Atlas provides five interactive modules:
| Module | Description |
|---|---|
| **UMAP Overview** | Visualize the global single-cell landscape colored by CellMajor / CellMin / Tissue / Cancer |
| **Gene Expression (Dotplot)** | Plot expression of any gene set across major cell types or sub-types |
| **Composition by Cancer** | Stacked bar chart showing cell-type proportions across cancer types |
| **Composition by Tissue** | Stacked bar chart showing cell-type proportions across tissues |
| **Gene in Subtypes** | Visualize gene expression across sub-clusters within a chosen major cell type |

## How to Use
This is a **standalone web app** — no installation, no server, no account needed.  
Just download it and open it on your own computer.
### Step 1. Download the repository
Click the green **<Code>** button at the top of this page → **Download ZIP**.  
Then unzip it anywhere on your computer (e.g. Desktop).
You will get a folder named `BoM-PanCancer-Atlas-main/` containing `index.html` and a `web_data/` folder.
### Step 2. Open it in your browser
>
> 1. Open a terminal / command prompt **inside the unzipped folder**.
> 2. Run:
>    ```bash
>    python -m http.server 8000
>    ```
>    (requires Python 3, which is pre-installed on Mac/Linux; on Windows, install from [python.org](https://www.python.org/))
> 3. Open your browser and go to **http://localhost:8000**
