# **UC Merced Woo Lab Dataset - Data Browser**  

This project is a **Shiny for Python** web application that visualizes **Zebrafish data** collected at **UC Merced's Woo Laboratory**. The application provides interactive plots, including:  

- **UMAP** (Uniform Manifold Approximation and Projection)  
- **Feature Plots**  
- **Dot Plots**  
- **Violin Plots**  

This tool is designed to help researchers explore and analyze single-cell data interactively.  

## **Installation**  

To run this application, you will need **Python 3.8 or later** installed on your system.  

### **1. Install Python (if not installed)**
- **Windows:** Download Python from [python.org](https://www.python.org/downloads/) and install it.  
- **Mac:** Python is usually pre-installed, but you can install/update it via Homebrew:  
  ```
  brew install python
  ```

### **2. Clone the Repository**
Run ```git clone https://github.com/your-repo-name/woo-lab-dataset-browser.git``` to download the project files.

### **3. Create a virtual environment (recommended)**
```
python3 -m venv venv  
source venv/bin/activate  # Mac/Linux  
venv\Scripts\activate     # Windows  

```

### **4. Install Dependencies**
```
pip install -r requirements.txt
```



## Running the Application
Once installed, you can start the application by running:
```
shiny run --host 0.0.0.0 --port 5000 app.py
```

**Access the app in your browser:**  
Open: `http://localhost:5000`

If you are running this on another computer in the same network, you can share the local IP:
```
http://your-ip-address:5000
```
To find your local IP:
- Windows: Run `ipconfig` and look for the "IPv4 Address."
- Mac/Linux: Run `ifconfig | grep "inet "` or `hostname -I`.

---

## Using the Data Browser
Once the application is running, you will see an interactive interface.

### Features:
1. **UMAP Plot:** View clusters of Zebrafish data points based on feature similarity.
2. **Feature Plots:** Display gene expression levels across different clusters.
3. **Dot Plots:** Compare expression levels across multiple conditions.
4. **Violin Plots:** Visualize distribution and variance of data features.


---

## Troubleshooting
### Common Issues & Fixes
| Issue | Possible Solution |
|--------|------------------|
| `ModuleNotFoundError: No module named 'shiny'` | Run `pip install shiny` |
| App fails to start due to missing dependencies | Run `pip install -r requirements.txt` |
| Unable to access from another computer | Use `--host 0.0.0.0`, check firewall settings |
| `Address already in use` error | Change the port: `shiny run --port 8080 app.py` |

---

## Credits
- **UC Merced Woo Laboratory**
- **Developers:** Brandon Jia
- **Special Thanks:** Woo Lab team for dataset collection

