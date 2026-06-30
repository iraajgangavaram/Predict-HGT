from pathlib import Path
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score


# ----------------------------
# DATA PATH
# ----------------------------

DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "processed" / "genome_ml_dataset.csv"


def main():

    print("\n=== Loading Dataset ===")
    print("Path:", DATA_PATH)

    df = pd.read_csv(DATA_PATH)

    print("\n=== Dataset Overview ===")
    print(df.shape)

    # ----------------------------
    # SAFETY CHECK (IMPORTANT)
    # ----------------------------

    if len(df) < 5:
        print("\n[ERROR] Dataset too small for ML training.")
        print("You need more genomes or samples.")
        return

    # ----------------------------
    # FEATURES / LABEL
    # ----------------------------

    if "label" not in df.columns:
        print("\n[ERROR] No 'label' column found.")
        return

    X = df.drop(columns=["label"])
    y = df["label"]

    # ----------------------------
    # TRAIN / TEST SPLIT SAFETY
    # ----------------------------

    if len(df) < 10:
        test_size = 0.2
    else:
        test_size = 0.25

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=42
    )

    # ----------------------------
    # MODEL
    # ----------------------------

    model = RandomForestClassifier(
        n_estimators=200,
        random_state=42
    )

    model.fit(X_train, y_train)

    # ----------------------------
    # PREDICTION
    # ----------------------------

    y_pred = model.predict(X_test)

    # ----------------------------
    # EVALUATION
    # ----------------------------

    print("\n=== Results ===")
    print("Accuracy:", accuracy_score(y_test, y_pred))

    print("\nClassification Report:\n")
    print(classification_report(y_test, y_pred, zero_division=0))


if __name__ == "__main__":
    main()