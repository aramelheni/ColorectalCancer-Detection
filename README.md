# 🧬 AI-Powered Gut Microbiome Analysis for Colorectal Cancer Detection  

Colorectal cancer (CRC) is one of the leading causes of cancer-related deaths worldwide.  
Early detection can significantly improve survival rates. Recent studies show that the **gut microbiome** plays a key role in CRC development, where microbial DNA patterns can serve as biomarkers for early diagnosis.  

This project provides an **end-to-end automated pipeline** for CRC detection from **gut microbiome DNA sequences**, consisting of two main phases:  

---

## 🔄 Project Phases  

### **Phase 1 Automated Data Processing Pipeline**  
- Traditionally, microbiome DNA preprocessing is performed with R scripts (e.g., DADA2), requiring manual steps.
 <img width="657" height="567" alt="image (3)" src="https://github.com/user-attachments/assets/6280be42-7003-435c-9816-3c56608f3ab4" />

- We automated this workflow entirely in **Python**, so the user only needs to provide the raw **DNA sequences (.fasta files)**.  
- The pipeline automatically:  
  - Processes and cleans the sequences  
  - Sorts the data into the required folder structure  
  - Prepares a filtered dataset ready for machine learning  

👉 **Result:** A fully automated and scalable preprocessing step, removing manual intervention.  
<img width="370" height="163" alt="image" src="https://github.com/user-attachments/assets/f1fe3715-8930-4882-a697-17e8a0963f03" />


---

### **Phase 2 Machine Learning for CRC Detection**  
- Using the preprocessed data, we apply **machine learning models** to classify samples as CRC-positive or CRC-negative.  
- **Steps:**  
  1. Extract k-mers from DNA sequences  
  2. Vectorize k-mers into numerical features with `CountVectorizer`  
  3. Split dataset into train/test (80/20)  
  4. Train multiple models: Logistic Regression, Random Forest, SVM, KNN  
  5. Evaluate models using Accuracy, Precision, Recall, F1-score, and Confusion Matrices  

👉 **Result:** Logistic Regression achieved the best performance, making it the most reliable model for deployment.  

---

## ⚙️ Tools & Libraries  

- **Python**: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`, `biopython`  
- **Automation**: R-to-Python pipeline integration  
- **Models**: Logistic Regression, Random Forest, SVM, KNN  

---

## 🖥️ Demo  

🎥 Watch the demo here: [YouTube Video](https://youtu.be/XXVqxhAETTk)  

---

**Supervised by:**  
Soumaya Jebara (UM6SS)  
Asma Amdouni (SMU)  

---

## ✅ Conclusion  

This project demonstrates that **automated DNA preprocessing + machine learning** can provide a scalable and reliable solution for early colorectal cancer detection from microbiome data.  
By automating the preprocessing pipeline and testing multiple models, we show that **Logistic Regression** is the most effective approach, paving the way for clinical integration.  

---

## License

This project is open-source and available under the MIT License. See the [LICENSE](LICENSE) file for details.
