def plot_confusion_matrix_custom(y_true, y_pred, title="Confusion Matrix", tag=""):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import confusion_matrix

    cm_ = confusion_matrix(y_true, y_pred)
    TN, FP, FN, TP = cm_.ravel()
    plot_cm = np.array([[TP, FP], [FN, TN]])

    cell_colors = np.array([
        ['seagreen', 'mediumseagreen'],
        ['darkseagreen', 'lightgreen']
    ])

    fig, ax = plt.subplots(figsize=(5, 4))
    for r in range(2):
        for c in range(2):
            ax.add_patch(plt.Rectangle((c, 1-r), 1, 1, color=cell_colors[r, c]))
            text_color = 'black' if (r, c) == (0, 0) else 'white'
            ax.text(c+0.5, 1.5-r, plot_cm[r, c],
                    ha='center', va='center', fontsize=12)

    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_xticks([0.5, 1.5])
    ax.set_xticklabels(['True Positive', 'True Negative'])
    ax.xaxis.tick_top()
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(['Predicted Negative', 'Predicted Positive'], rotation=90)

    ax.set_title(title)

    # filename based on tag
    if tag == "all_benchmark":
        filename = "CM_benchmark_ALL_features.png"
    else:
        filename = "CM_benchmark_selected.png"

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()
    print(f"Saved confusion matrix to: {filename}")
