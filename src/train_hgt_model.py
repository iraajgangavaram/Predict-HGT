from pathlib import Path
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score


def main():

    print("\n=== Training HGT Predictor ===\n")

    df = pd.read_csv("data/processed/hgt_dataset.csv")

    print("Dataset shape:", df.shape)

    # -------------------------
    # FEATURES / LABEL
    # -------------------------

    X = df[["similarity"]]
    y = df["label"]

    # -------------------------
    # TRAIN TEST SPLIT
    # -------------------------

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.25,
        random_state=42
    )

    # -------------------------
    # MODEL
    # -------------------------

    model = RandomForestClassifier(
        n_estimators=200,
        random_state=42
    )

    model.fit(X_train, y_train)

    # -------------------------
    # PREDICTIONS
    # -------------------------

    y_pred = model.predict(X_test)

    # -------------------------
    # RESULTS
    # -------------------------

    print("\n=== Results ===")
    print("Accuracy:", accuracy_score(y_test, y_pred))

    print("\nClassification Report:\n")
    print(classification_report(y_test, y_pred, zero_division=0))


if __name__ == "__main__":
    main()