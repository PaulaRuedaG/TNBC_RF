# -*- coding: utf-8 -*-
"""RF_Model.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1zF4LOXf8NV5EQl3xayV-ZUhLQ1hhy9zD
"""

#Import libraries
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix, roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.model_selection import cross_val_score

#Import data
datos_ml = pd.read_csv("datos_ML_Conjuntos_1.csv",index_col=0)
print(datos_ml.shape)
X = datos_ml.drop('grupos', axis=1)
y = datos_ml['grupos']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
print(X_train.shape, X_test.shape)

#Create Random Forest Model
modelo_rf = RandomForestClassifier(n_estimators=200, random_state=42, class_weight='balanced')
modelo_rf.fit(X_train, y_train)

y_pred = modelo_rf.predict(X_test)

print("Confussion matrix:")
print(confusion_matrix(y_test, y_pred))

print("\nClassification report:")
print(classification_report(y_test, y_pred))

print("\nModel precision:")
print(accuracy_score(y_test, y_pred))

pos_label = 'pCR'  # Specify the positive class label

# Get the probabilities of the positive class
y_probs = modelo_rf.predict_proba(X_test)[:, 1]

# Calculate the ROC curve
fpr, tpr, _ = roc_curve(y_test, y_probs, pos_label=pos_label)  # Pass pos_label
roc_auc = auc(fpr, tpr)

# Plot the ROC curve
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label=f'AUROC = {roc_auc:.2f}')
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.title('ROC Curve')
plt.legend(loc="lower right")
plt.show()

print(f"Área bajo la curva ROC (AUROC): {roc_auc:.4f}")
plt.savefig("AUROC.jpg", dpi=300)  # para un journal es mejor TIFF a 300dpi
plt.show()

# Calcular curva PR
# Pass the pos_label to precision_recall_curve
precision, recall, _ = precision_recall_curve(y_test, y_probs, pos_label='pCR')
auprc = average_precision_score(y_test, y_probs, pos_label='pCR')

# Graficar la curva PR
plt.figure()
plt.plot(recall, precision, color='green', lw=2, label=f'AUPRC = {auprc:.2f}')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc="upper right")
plt.show()

print(f"Área bajo la curva PR (AUPRC): {auprc:.4f}")
plt.savefig("AUPRC.jpg", dpi=300)  # para un journal es mejor TIFF a 300dpi
plt.show()

modelo_rf = RandomForestClassifier(n_estimators=200, random_state=42, class_weight='balanced')
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
cv_scores = cross_val_score(modelo_rf, X, y, cv=cv)

#Cross validation
cv_scores = cross_val_score(modelo_rf, X, y, cv=10)
print("Puntuaciones de validación cruzada en cada fold:", cv_scores)
print("Puntuación media:", cv_scores.mean())