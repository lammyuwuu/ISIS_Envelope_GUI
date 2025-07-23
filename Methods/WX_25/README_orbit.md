v   st# 📦 BPM Simulation for Accelerator Physics

This guide explains how to install Python dependencies and run the Streamlit app from `app.py`.

---

## 🔧 Prerequisites

Make sure you have the following installed:

- Python **3.7 or higher**
- `pip` (Python package installer)
- (Optional) A virtual environment tool like `venv` or `conda`

---

## 📥 Step 1: Clone the Repository

```bash
git clone git@https://gitlab.stfc.ac.uk/vfl94988/synchrotron_model
cd .\synchrotron_model
```
## OR
Download from gitlab repository using Code drop down arrow, then Download Directory and select zip.
Then in Downloads Extract all to a file in your Documents folder
---

## 📦 Step 2: Install Dependencies

Run this in your project root where `requirements.txt` is located:

```bash
pip install -r requirements.txt
```

### 💡 Using a virtual environment (recommended):

```bash
# Create a virtual environment
python -m venv venv

# Activate it
source venv\Scripts\activate        # On Linux: venv/bin/activate  

# Install dependencies
pip install -r requirements.txt
```

---
## Step 3: Run Jupyter Notebooks

Run the Jupyter Notebook 00_Plot_Beam_Orbit.ipynb
Then run the Jupyter Notebook 01_Correct_Orbit.ipynb

## 🚀 Step 4: Run the Streamlit App

To launch the app:

```bash
cd .\ISIS_Synchrotron_Model\Methods\WX_25\
streamlit run app.py
```

This will open the app in your browser at [http://localhost:8501](http://localhost:8501).

---

## 🛠 Common Issues

- **`ModuleNotFoundError`**: Double-check that dependencies were installed correctly.
- **Streamlit doesn’t open in browser**: Manually open [http://localhost:8501](http://localhost:8501).

---

## 🧹 Optional: Deactivate Virtual Environment
First stop streamlit
Use ctrl+c to stop, then
```bash
deactivate
```