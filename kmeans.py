import numpy as np

def kmeans(points, k, max_iter=300, eps=0.001):
    n = len(points)
    d = len(points[0])
    
    # Initialize centroids
    centroids = points[:k]
    
    clusters = [[] for _ in range(k)]
    
    def dist(p1, p2):
        return np.linalg.norm(np.array(p1) - np.array(p2))
    
    def assign_clusters():
        for i in range(k):
            clusters[i] = []
        
        for point in points:
            nearest = float('inf')
            idx = 0
            c_idx = 0
            for centroid in centroids:
                dis = dist(point, centroid)
                if dis < nearest:
                    nearest = dis
                    c_idx = idx
                idx += 1
            clusters[c_idx].append(point)
    
    def update_centroids():
        nonlocal centroids
        flag = 0
        new_centroids = []
        
        for cluster in clusters:
            if len(cluster) == 0:
                new_centroids.append([0.0] * d)
                continue
            new_centroid = np.mean(np.array(cluster), axis=0).tolist()
            new_centroids.append(new_centroid)
        
        for i in range(k):
            if dist(centroids[i], new_centroids[i]) >= eps:
                flag = 1
        
        centroids = new_centroids
        return flag
    
    for _ in range(max_iter):
        assign_clusters()
        if update_centroids() == 0:
            break
    
    # Assign labels
    labels = np.zeros(n, dtype=int)
    for i, point in enumerate(points):
        nearest = float('inf')
        idx = 0
        c_idx = 0
        for centroid in centroids:
            dis = dist(point, centroid)
            if dis < nearest:
                nearest = dis
                c_idx = idx
            idx += 1
        labels[i] = c_idx
    
    return labels.tolist()