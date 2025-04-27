import sys
    
if len(sys.argv) not in [3,4]:
    print("An Error Has Occurred")
    sys.exit()

if len(sys.argv)==3:
    iter = 200

    f = sys.argv[2]
    f2 = f[-4:]
    if f2 != ".txt":
        print("NA")
        sys.exit() 

else:
    iter = sys.argv[2]

    f = sys.argv[3]
    f2 = f[-4:]
    if f2 != ".txt":
        print("NA")
        sys.exit() 

    if iter.isdigit() == False or int(iter) <= 1 or int(iter) >= 1000 or iter[0]=="0":
        print ("Invalid maximum iteration!")
        sys.exit() 

f1 = open(f)

N = 0
for line in f1:
    N += 1

f1 = open(f)
d = len(f1.readline().split(","))


eps = 0.001


k = sys.argv[1]

if k.isdigit() == False or int(k) <= 1 or int(k) >= N or k[0]=="0":
    print ("Invalid number of clusters!")
    sys.exit() 


f1 = open(f)
points = []
for line in f1:
    points.append(line.split(","))
   
centroids = points[0:int(k)]

clusters = []
for i in range(int(k)):
    clusters.append([])


def dist(point1 , point2):
    sum = 0
    for i in range(d):
        sum += (float(point1[i])-float(point2[i]))**2
    return sum**0.5
    


def point_to_cluster():
    for i in range(len(clusters)):
        clusters[i] = []

    for point in points:
        nearest = float('inf')
        idx = 0
        c_idx = 0
        for centroid in centroids:
            dis = dist(point , centroid)
            if dis < nearest:
                nearest = dis
                c_idx = idx
            idx += 1

        clusters[c_idx].append(point)


def cluster_avg():
    idx = 0
    flag = 0
    for cluster in clusters:
        tmp = []
        l = len(cluster)
        
        if l==0:
            continue
 
        for j in range(d):
            sum = 0
            i = 0
            while i < l:
                sum += float(cluster[i][j])
                i += 1
            sum = sum / l
            tmp.append(sum)
            j += 1
        
        tmp1 = centroids[idx]
        centroids[idx] = tmp
        
        idx += 1

        if len(tmp) != 0 and len(tmp1) !=0:
            if dist(tmp1, tmp) > eps:
                flag = 1
    
            
    
    if flag == 0:
        return False


for i in range(int(iter)):
    try:
        point_to_cluster()
        if cluster_avg() == False:
            break
    except:
        print("An Error Has Occurred")
        sys.exit() 

    
    
for centroid in centroids:
    print(",".join(f"{x:.4f}" for x in centroid))