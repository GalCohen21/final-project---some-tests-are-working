import sys
import numpy as np
import symnmf

def main():
    # command line arguments check
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)

    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    valid_goals = ['symnmf', 'sym', 'ddg', 'norm']
    if goal not in valid_goals:
        print("An Error Has Occurred")
        sys.exit(1)

    # read file
    try:
        X = np.loadtxt(file_name, delimiter=',')
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    
    X_list = X.tolist()
    
    # Calculate full symNMF
    if goal == "symnmf":
        # Compute W
        W = symnmf.norm(X_list)
        
        np.random.seed(1234)
        # Initialize H
        n = len(X)
        m = np.mean(W)
        H = np.random.uniform(0, 2 * np.sqrt(m/k), (n, k)).tolist()
        
        # Calculate H
        result = symnmf.symnmf(H, W)
            
    # Calculate similarity matrix
    elif goal == "sym":
        result = symnmf.sym(X_list)
    
    # Calculate diagonal degree matrix
    elif goal == "ddg":
        result = symnmf.ddg(X_list)
    
    # Calculate normalized similarity matrix
    elif goal == "norm":
        result = symnmf.norm(X_list)

    for row in result:
        print(','.join(f"{x:.4f}" for x in row))

if __name__ == "__main__":
    main()