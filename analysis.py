import sys
import numpy as np
from sklearn.metrics import silhouette_score
import symnmf
import kmeans

def main():
    if len(sys.argv) != 3:
        print("An Error Has Occurred")
        sys.exit(1)
    
    try:
        k = int(sys.argv[1])
        file_name = sys.argv[2]
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    
    try:
        data = np.loadtxt(file_name, delimiter=',')
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    
    n = data.shape[0]
    
    if k < 1 or k >= n:
        print("An Error Has Occurred")
        sys.exit(1)
    
    data_list = data.tolist()
    
    W = symnmf.norm(data_list)
    
    np.random.seed(1234)
    m = np.mean(W)
    H = np.random.uniform(0, 2 * np.sqrt(m/k), (n, k)).tolist()
    
    H_final = symnmf.symnmf(H, W)
    H_np = np.array(H_final)
    
    nmf_labels = np.argmax(H_np, axis=1)
    
    # --- Use your own KMeans! ---
    kmeans_labels = kmeans.kmeans(data_list, k, max_iter=300)
    
    try:
        nmf_score = silhouette_score(data, nmf_labels)
        kmeans_score = silhouette_score(data, kmeans_labels)
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    
    print(f"nmf: {nmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()