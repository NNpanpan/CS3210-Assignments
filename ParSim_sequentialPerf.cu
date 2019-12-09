#include <stdio.h>
#include <string.h>
#include <bits/stdc++.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace std;

const double eps = 1e-8;
const double PI = acos(-1);

__constant__ double EPS;

__device__ long cnt;
int host_n;
int step, host_s;
double host_l, host_r;


/*
-----STRUCTURES-----
*/
struct Point
{ 
    double x, y;

    __host__ __device__ Point(double cx, double cy) {
        x = cx;
        y = cy;
    }

    __host__ __device__ Point operator + (const Point &oth) {
        return Point(x + oth.x, y + oth.y);
    }
    __host__ __device__ Point operator - (const Point &oth) {
        return Point(x - oth.x, y - oth.y);
    }
    __host__ __device__ Point operator * (double k) {
        return Point(x*k, y*k);
    }

    __host__ __device__ double operator * (const Point &oth) {
        return x * oth.x + y * oth.y;
    }
};

typedef Point Velocity;

struct Particle 
{
    struct Point center = Point(0.0, 0.0);
    double radius;
    Velocity v = Velocity(0.0, 0.0);
    long id;

    long wColl;
    long pColl;

    double lastCol;

    __host__ __device__ Particle(struct Point C, double R, Velocity V, long ID) {
        center = C;
        radius = R;
        v = V;
        id = ID;
        wColl = 0;
        pColl = 0;
        lastCol = -1.0;
    }

    __host__ __device__ void Move(double mtime, double limit){
        double EPS = 1e-8;
        double newX = center.x + v.x * mtime;
        if (newX >= limit - radius + EPS) newX = limit-radius;
        else if (newX + EPS <= radius) newX = radius;
        
        double newY = center.y + v.y * mtime;
        if (newY >= limit - radius + EPS) newY = limit-radius;
        else if (newY + EPS <= radius) newY = radius;

        center = Point(newX, newY);
    } //only check for wall "touching"

    __host__ __device__ bool Overlap(const Particle &oth) {
        double EPS = 1e-8;
        Point shift = center - oth.center;
        double sqDist = shift * shift;
        double sqSumR = (radius + oth.radius) * (radius + oth.radius);
        return sqDist <= sqSumR - EPS; 
    }
};

struct Collision 
{
    long par1, par2;
    double ctime;
    long checkSum;

    __host__ __device__ Collision(long p1, long p2, double t, long sum = 0) {
        par1 = p1 > p2 ? p1 : p2;
        par2 = p1 + p2 - par1;
        ctime = t;
        checkSum = sum;
    }

    __host__ __device__ bool operator < (const Collision &oth) const {
        if (ctime > oth.ctime + eps)
            return false;
        else if (ctime + eps < oth.ctime)
            return true;
        else {
            if (par2 == oth.par2) return par1 < oth.par1;
            return par2 < oth.par2;
        }
    }
};

__managed__ Collision* colls;

struct collHeap 
{
    long n;

    __host__ collHeap() {
        n = 0;
    }

    __host__ __device__ void shiftUp(long index) {
        long k = index;
        while (k > 1 && colls[k] < colls[k/2]) {
            Collision temp = colls[k];
            colls[k] = colls[k/2];
            colls[k/2] = temp;
            k /= 2;
        }
    }

    __host__ __device__ void shiftDown(long index) {
        long k = index;
        while (k*2 <= n) {
            long dir = k * 2;
            if (k * 2 + 1 <= n && colls[k * 2 + 1] < colls[k * 2])
                dir = k * 2 + 1;
            if (colls[k] < colls[dir]) break;
            Collision temp = colls[k];
            colls[k] = colls[dir];
            colls[dir] = temp;
            k = dir;
        }
    }

    __host__ void insert(Collision C) {
        n++;
        colls[n] = C;
        shiftUp(n);
    }

    __host__ __device__ Collision pop() {
        Collision ret = colls[1];
        colls[1] = colls[n];
        n--;
        shiftDown(1);
        return ret;
    }

    __host__ bool isEmpty() {
        return n == 0;
    }

    __host__ void build() {
        for (long k = n/2; k >= 1; k--) {
            shiftDown(k);
        }
    }


	__host__ void update(long newPos) {
		for (long k = n; k <= newPos; k++) {
			shiftUp(k);
		}
	}
};


/*
-----COLLISION TIME-----
*/
__host__ __device__ double CollideTimePar(Particle &a, Particle &b) {
    double EPS = 1e-8;
	Velocity diffV = a.v - b.v;
    struct Point diffC = a.center - b.center;
    double distSq = diffC * diffC;
    double touchD = 4*a.radius*a.radius;
    if (distSq + EPS < touchD) {
		return 0.0;
    }

    double bVal = diffV * diffC;
    double aVal = diffV * diffV;
    double delta = bVal*bVal - aVal*(distSq - touchD);
    if (delta <= EPS) {
        return -2.0;
    } else {
        double ret1 = (-bVal - sqrt(delta))/aVal;
        double ret2 = (-bVal + sqrt(delta))/aVal;
        return ret1 >= 0.0 + EPS ? ret1 : ret2;
    }
}

__host__ __device__ double CollideTimeParPerf(Particle &a, Particle &b, double cur) {
    double EPS = 1e-8;
	Velocity diffV = a.v - b.v;
    struct Point diffC = a.center - b.center;
    double distSq = diffC * diffC;
    double touchD = 4*a.radius*a.radius;
    if (distSq + EPS < touchD) {
		return -2.0;
    }

    double bVal = diffV * diffC;
    double aVal = diffV * diffV;
    double delta = bVal*bVal - aVal*(distSq - touchD);
    if (delta <= EPS) {
        return -2.0;
    } else {
        double ret1 = (-bVal - sqrt(delta))/aVal + cur;
        double ret2 = (-bVal + sqrt(delta))/aVal + cur;
        double maxLastCol = max (a.lastCol, b.lastCol);
        if (ret1 + EPS < maxLastCol && maxLastCol + EPS < ret2) return maxLastCol;
        else if (maxLastCol + EPS < ret1) return ret1;
        else return -2.0;
    }
}

__host__ __device__ double CollideTimeWall(Particle &a, double boxLength) {
    double r = a.radius;

    double EPS = 1e-8;
    double xTime, yTime;
    if (a.v.x > 0.0 + EPS) {
        xTime = abs(boxLength - r - a.center.x)/a.v.x;
    } else if (a.v.x <= -EPS) {
        xTime = abs((a.center.x - r)/a.v.x);
    } else {
        xTime = DBL_MAX;
    }

    if (a.v.y > 0.0 + EPS) {
        yTime = abs(boxLength - r - a.center.y)/a.v.y;
    } else if (a.v.y <= -EPS) {
        yTime = abs((a.center.y - r)/a.v.y);
    } else {
        yTime = DBL_MAX;
    }

    return xTime > yTime + EPS ? yTime : xTime;
}



/*
-----CONDITIONS TESTING-----
*/
__host__ bool isValidCollisionTime (double ctime, double curTime, double lim) 
{
    return ctime + eps > curTime && ctime + eps < lim;
}

__host__ bool isValidChecksum (Particle* host_pars, long id1, long id2, long checkSum) 
{
    if (id2 == -1) {
        return host_pars[id1].wColl + host_pars[id1].pColl == checkSum;
    } else {
        return host_pars[id1].wColl + host_pars[id1].pColl
            + host_pars[id2].wColl + host_pars[id2].pColl == checkSum;
    }
}

__host__ bool validLastCols (Particle* pars, long id1, long id2, double ctime) 
{
    if (id2 == -1) {
        return pars[id1].lastCol <= ctime + eps;
    } else {
        return pars[id1].lastCol <= ctime + eps && pars[id2].lastCol <= ctime + eps;
    }
}

__host__ bool wallColX (int ids, Particle* pars) 
{
    return abs(pars[ids].center.x - pars[ids].radius) <= eps 
        || abs(host_l - pars[ids].center.x - pars[ids].radius) <= eps;
}

__host__ bool wallColY (int ids, Particle* pars) 
{
    return abs(pars[ids].center.y - pars[ids].radius) <= eps 
        || abs(host_l - pars[ids].center.y - pars[ids].radius) <= eps;
}



/*
-----COLLISION RESOLVING-----
*/
__host__ void parCol(Particle* host_pars, int a, int b) 
{
    struct Point un = host_pars[a].center - host_pars[b].center;
    double unLen = sqrt(un * un);
    un = Point (un.x / unLen, un.y / unLen);
    struct Point ut = Point(-un.y, un.x);

    double v1n = un * host_pars[a].v;
    double v1t = ut * host_pars[a].v;
    double v2n = un * host_pars[b].v;
    double v2t = ut * host_pars[b].v;

    // db v1nf = v2n; db v2nf = v1n;
    host_pars[a].v = Point(v2n * un.x + v1t * ut.x, v2n * un.y + v1t * ut.y);

    host_pars[b].v = Point(v1n * un.x + v2t * ut.x, v1n * un.y + v2t * ut.y);

    host_pars[a].pColl++;
    host_pars[b].pColl++;
}

__host__ void wallCol(Particle* host_pars, int p) 
{
    if (wallColX(p, host_pars)) 
        host_pars[p].v = Velocity(-host_pars[p].v.x, host_pars[p].v.y);
        
    if (wallColY(p, host_pars)) 
        host_pars[p].v = Velocity(host_pars[p].v.x, -host_pars[p].v.y);

    host_pars[p].wColl++;
}

long totalP = 0; // to count total collisions between particles
long totalW = 0; // to count total wall collisions
/*
-----PRINT MODE-----
*/
__global__ void simulate_step(long n, double l, Particle* device_pars, long lim)
{
	long len = n;
    long i = blockIdx.x * lim;
	long k = threadIdx.x * lim;
    
	//printf("kernel on thread: %ld\n", i+k);
    for (long cnt1 = 0; i + cnt1 < len && cnt1 < lim; cnt1 ++) {
        for (long cnt2 = 0; k + cnt2 < len && cnt2 < lim; cnt2 ++) {
			// if (i+cnt1 < k+cnt2) printf("particles %ld and %ld\n", i+cnt1, k+cnt2);
            if (i + cnt1 == k + cnt2) {
                double wallTime = CollideTimeWall(device_pars[i + cnt1], l);
                
				if (wallTime + EPS < 0 || wallTime > 1.0 + EPS)
                    continue;
                // printf("%ld will hit wall at %10.8f\n", i+cnt1, wallTime);
                long j = atomicAdd((unsigned long long*) &cnt, (unsigned long long) 1);
                // printf("collision %ld\n", j);
                colls[j] = Collision( i + cnt1, -1, wallTime );
            } else {
                long p1 = i + cnt1;
                long p2 = k + cnt2;
                double colTime = CollideTimePar(device_pars[p1], device_pars[p2]);
                if (colTime + EPS < 0 || colTime > 1.0 + EPS) {
                    continue;
                }
                
                long j = atomicAdd((unsigned long long*) &cnt, (unsigned long long) 1);
                // printf("collision %ld\n", j);
                colls[j] = Collision( p1, p2, colTime );
            }
        }
    }
    
}

__host__ void step_resolve(Particle* host_pars, Particle* device_pars, long coll_size) 
{
    // cout << "sorting size " << coll_size << "\n";
    std::sort(colls, colls + coll_size);
    vector<int> hasCollision (host_n, 0);
	
	// cout << "resolving\n";
    for(long i = 0; i < coll_size; i++) {
        Collision c = colls[i];
        if (c.par2 == -1 && !hasCollision[c.par1]) {
            ///wall collision

			long ID = c.par1;
			// printf("particle %ld hit wall at %10.8f\n", ID, c.ctime);
            host_pars[ID].Move(c.ctime, host_l);
            wallCol(host_pars, ID);

            hasCollision[ID] = 1;

            if (1.0 - c.ctime > eps)
                host_pars[ID].Move(1.0-c.ctime, host_l);

            totalW++;
        }
        else if (!hasCollision[c.par1] && !hasCollision[c.par2]) {
            ///2 particles

            long a = c.par1;
            long b = c.par2;
			// printf("particles %ld and %ld hit at %10.8f\n", a, b, c.ctime);
            host_pars[a].Move(c.ctime, host_l);
            host_pars[b].Move(c.ctime, host_l);
            parCol(host_pars, a, b);

            hasCollision[a] = 1;
            hasCollision[b] = 1;

            if (1.0 - c.ctime > eps) {
                host_pars[a].Move(1.0-c.ctime, host_l);
                host_pars[b].Move(1.0-c.ctime, host_l);
            }

            totalP++;
        }
    }
    // cout << "done resolving\n";
    for(long i = 0; i < host_n; i++) if (!hasCollision[i]) {
        host_pars[i].Move(1.0, host_l);
    }
    // cout << "done moving. new loop ready\n";
}

/*
-----PERF MODE-----
*/
__global__ void simulate_perf(long n, long l, Particle* device_pars, long lim, double s) 
{
    long len = n;
    long i = blockIdx.x * lim;
    long k = threadIdx.x * lim;
      
    for (long cnt1 = 0; i + cnt1 < len && cnt1 < lim; cnt1 ++) {
        for (long cnt2 = 0; k + cnt2 < len && cnt2 < lim; cnt2 ++) {
            // if (i+cnt1 < k+cnt2) printf("particles %ld and %ld\n", i+cnt1, k+cnt2);
            if (i + cnt1 == k + cnt2) {
                long p = i+cnt1;
                double wallTime = CollideTimeWall(device_pars[p], l);
                
				if (wallTime + EPS < 0 || wallTime > s + EPS)
                    continue;
                
                long pos = atomicAdd((unsigned long long*) &cnt, (unsigned long long) 1)+1;
                // printf("collision %ld\n", pos);
                colls[pos] = Collision( p, -1, wallTime );
            } else {
                long p1 = i + cnt1;
                long p2 = k + cnt2;
                if (p1 > p2) continue;
                // printf("particles %ld and %ld\n", p1, p2);

                double colTime = CollideTimeParPerf(device_pars[p1], device_pars[p2], 0.0);
                if (colTime + EPS < 0 || colTime > s + EPS) {
                    continue;
                }
                
                long pos = atomicAdd((unsigned long long*) &cnt, (unsigned long long) 1)+1;
                // printf("collision %ld\n", pos);
                colls[pos] = Collision( p1, p2, colTime );
            }
        }
    }
}

__host__ void hostRecalcSingle(Particle* pars, long id, double step, long n, long l, double cur, collHeap* heap)
{
    long i = 0;
    double EPS = 1e-8;
    for (long ct = 0; i + ct < n; ct++) {
        if (i + ct == id) {
            double wallTime = CollideTimeWall(pars[id], l) + cur;
                
            if (wallTime + EPS < cur || wallTime > step + EPS)
                continue;
            
            if (wallTime <= pars[id].lastCol + EPS)
                wallTime = pars[id].lastCol;

            long sum = pars[id].wColl + pars[id].pColl;
            heap->insert(Collision( id, -1, wallTime, sum ));
        } else {
            double colTime = CollideTimeParPerf(pars[i+ct], pars[id], cur);
            if (colTime + EPS < cur || colTime > step + EPS) {
                continue;
            }
            long sum = pars[id].wColl + pars[id].pColl
                    + pars[i+ct].wColl + pars[i+ct].pColl;
            heap->insert(Collision( id, i+ct, colTime, sum ));
        }
    }
}

__host__ void hostRecalcDouble(Particle* pars, long id1, long id2, double step, long n, long l, double cur, collHeap *heap)
{
    long i = 0;
    double EPS = 1e-8;
    for (long ct = 0; i + ct < n; ct++) {
        long p = i + ct;
        if (p == id1) {
            double wallTime = CollideTimeWall(pars[id1], l) + cur;
                
            if (wallTime + EPS < cur || wallTime > step + EPS)
                continue;
            
            if (wallTime <= pars[id1].lastCol + EPS)
                wallTime = pars[id1].lastCol;

            long sum = pars[id1].wColl + pars[id1].pColl;
            heap->insert(Collision( id1, -1, wallTime, sum ));
        } else if (p == id2) {
            double wallTime = CollideTimeWall(pars[id2], l) + cur;
                
            if (wallTime + EPS < cur || wallTime > step + EPS)
                continue;

            if (wallTime <= pars[id2].lastCol + EPS)
                wallTime = pars[id2].lastCol;
            
            long sum = pars[id2].wColl + pars[id2].pColl;
            heap->insert(Collision( id2, -1, wallTime, sum ));
        } else {    
            double colTime1 = CollideTimeParPerf(pars[p], pars[id1], cur);
            double colTime2 = CollideTimeParPerf(pars[p], pars[id2], cur);
            
            if (colTime1 + EPS >= cur && colTime1 <= step + EPS) {
                long sum = pars[id1].wColl + pars[id1].pColl
                    + pars[p].wColl + pars[p].pColl;
                heap->insert(Collision( id1, p, colTime1, sum ));
            }

            if (colTime2 + EPS >= cur && colTime2 <= step + EPS) {
                long sum = pars[id2].wColl + pars[id2].pColl
                    + pars[p].wColl + pars[p].pColl;
		        heap->insert(Collision( id2, p, colTime2, sum ));
            }
            
        }
    }
}

__host__ void perf_resolve(Particle* host_pars, Particle* device_pars_ptr, double step) 
{
    thrust::device_ptr<Particle> device_pars(device_pars_ptr);
    double curTime = 0.0;
    double lapse;
    collHeap heap = collHeap();
    cudaMemcpyFromSymbol(&heap.n, cnt, sizeof(long));
    heap.build();
    // cout << "heap size: " << heap.n << "\n";
    long k = host_n >= 1e6 ? 4096 : host_n >= 1e5 ? 2048 : host_n >= 1e4 ? 512 : 256;
    long split = (host_n + k-1)/k;

    while (!heap.isEmpty()) {
        Collision col = heap.pop();
        if (!isValidCollisionTime(col.ctime, curTime, step) 
            || !isValidChecksum(host_pars, col.par1, col.par2, col.checkSum)
            || !validLastCols(host_pars, col.par1, col.par2, col.ctime))
            continue;

        lapse = col.ctime - curTime;
        curTime = col.ctime;
        //cudaMemcpyToSymbol(cnt, &heap.n, sizeof(long));
        // printf("collision at %10.8f\n", col.ctime);

        for (long i = 0; i < host_n; i++) {
            host_pars[i].Move(lapse, host_l);
        }

        if (col.par2 == -1) {
            long id = col.par1;
            wallCol(host_pars, id);
            host_pars[id].lastCol = ceil(col.ctime + eps);
            totalW++;

            //thrust::copy(host_pars, host_pars + host_n, device_pars);
            
            // re-calc
            
            hostRecalcSingle(host_pars, id, step, host_n, host_l, curTime, &heap);
        } else {
            long id1 = col.par1;
            long id2 = col.par2;
            parCol(host_pars, id1, id2);
            host_pars[id1].lastCol = ceil(col.ctime + eps);
            host_pars[id2].lastCol = ceil(col.ctime + eps);
            totalP++;

            //thrust::copy(host_pars, host_pars + host_n, device_pars);

            // re-calc

            hostRecalcDouble(host_pars, id1, id2, step, host_n, host_l, curTime, &heap);

            double ct = CollideTimePar(host_pars[id1], host_pars[id2]) + curTime;
            if (ct + eps >= 0 && ct <= step + eps) {
                long sum = host_pars[id1].wColl + host_pars[id1].pColl
                        + host_pars[id2].wColl + host_pars[id2].pColl;
                heap.insert(Collision( ct, id1, id2, sum ));
            }
        }
    }

    // printf("all collisions done at %10.8f, remaining %10.8f\n", curTime, step-curTime);

    if (step - curTime > eps) {
        // printf("moving everything for %10.8f\n", step-curTime);
        for (long i = 0; i < host_n; i++) {
            host_pars[i].Move(step - curTime, host_l);
        }
    }
}



/*
-----STATISTICS-----
*/
__host__ void print_particles(int step, Particle* host_pars)
{
    int i;
    for (i = 0; i < host_n; i++) {
        Particle A = host_pars[i];
        printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", step, i, A.center.x, A.center.y,
            A.v.x, A.v.y);
    }
}

__host__ void print_statistics(int num_step, Particle* pars)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f %ld %ld\n", num_step, i, pars[i].center.x,
            pars[i].center.y, pars[i].v.x, pars[i].v.y,
            pars[i].pColl, pars[i].wColl);
    }
}


/*
-----RANDOMIZER-----
*/
__host__ double rng01(mt19937 &myRng) 
{
    unsigned value = myRng();
    return (value + .0) / myRng.max();
}

__host__ Particle RandomParticle(double boxLength, double radius, int index, 
    mt19937& myRng) 
    {
    
    double px = radius + (boxLength - 2 * radius) * rng01(myRng);
    double py = radius + (boxLength - 2 * radius) * rng01(myRng);

    double speedMin = boxLength / (8 * radius); 
    double speedMax = boxLength / 4;
    assert(speedMin <= speedMax);

    double speed = speedMin + (speedMax - speedMin) * rng01(myRng);
    double speedAng = rng01(myRng) * (2 * PI);

    double vx = cos(speedAng) * speed; 
    double vy = sin(speedAng) * speed;
    return Particle(Point(px, py), radius, Velocity(vx, vy), index); 
}

__host__ bool AddParticle(const Particle &part, int pos, Particle* pars)
{
	// cout << "checking particle " << pos << "\n";
    for(int i = 0; i < pos; i++) if (pars[i].Overlap(part)) return 0; 

   	pars[pos] = part;
 
    return 1;
}


int main(int argc, char** argv)
{
    time_t start, end;

    string mode;

    cin >> host_n >> host_l >> host_r >> host_s >> mode;

    time(&start);
    
    long K = host_n <= 1e5 + 50000 ? 1 << max(5, (int)floor(log2(host_n)) / 2) : 1 << 10;
    
    long num_blocks = (host_n + K-1)/K;
    long num_threads = num_blocks;

    thrust::host_vector<Particle> host_pars(host_n, Particle(Point(-1,-1), host_r, Velocity(-1,-1),-1));
    thrust::device_vector<Particle> device_pars(host_n, Particle(Point(-1,-1), host_r, Velocity(-1,-1),-1));  
    cudaMallocManaged((void**)&colls, sizeof(Collision) * host_n * host_n);
    Particle* host_pars_ptr = thrust::raw_pointer_cast(&host_pars[0]);
    Particle* device_pars_ptr = thrust::raw_pointer_cast(&device_pars[0]);


    long idx;
    double x, y, vx, vy;
    if (scanf("%ld %lf %lf %lf %lf", &idx, &x, &y, &vx, &vy) != EOF) {
        // cout << "read from input\n";
        long i = 0;
        host_pars[i] = (Particle(Point(x,y),
                    host_r,
                    Velocity(vx, vy),
                    idx));
        for(i = 1; i < host_n; i++) {
            
            cin >> idx >> x >> y >> vx >> vy;
            host_pars[i] = (Particle(Point(x,y),
                    host_r,
                    Velocity(vx, vy),
                    idx));
        }
    } else {
		// cout << "randomize\n";
        string seed = "agony";
        if (argc == 2) {
            seed = string(argv[1]);
            seed_seq ss (seed.begin(), seed.end());
            mt19937 myRng(ss);
            for(int i = 0; i < host_n; i++) {
                Particle newP = RandomParticle(host_l, host_r, i, myRng); 
                if (!AddParticle(newP, i, host_pars_ptr)) i--;
            }
        } else {
            mt19937 myRng(time(NULL));

            for(int i = 0; i < host_n; i++) {
                Particle newP = RandomParticle(host_l, host_r, i, myRng); 
                if (!AddParticle(newP, i, host_pars_ptr)) i--;
            }
        }   
    }

    /* Copy to GPU constant memory */
    cudaMemcpyToSymbol(EPS, &eps, sizeof(eps));


    long total_colls;
    
    print_particles(0, host_pars_ptr);

    if (mode == "print") {
        // cout << "print mode\n";
        
        for (step = 0; step < host_s; step++) {
			// cout << "step: " << step << "\n";
            if (step > 0) print_particles(step, host_pars_ptr);
           
            // device_pars = host_pars;
            thrust::copy(host_pars.begin(), host_pars.end(), device_pars.begin());
            // cout << "copied to device_pars\n";

            total_colls = 0;
            cudaMemcpyToSymbol(cnt, &total_colls, sizeof(total_colls));
            
            // cout << "set cnt to 0\n";

            /* Call the kernel */
            // cout << "got pointers, calling kernel\n";
            simulate_step<<<num_blocks, num_threads>>>(host_n, host_l, device_pars_ptr, K);

            /* Barrier */
            cudaDeviceSynchronize();
			// cout << "done calculation\n";           	
 
            cudaMemcpyFromSymbol(&total_colls, cnt, sizeof(total_colls));
            // cout << "copied cnt back to total_colls\n";

            step_resolve(host_pars_ptr, device_pars_ptr, total_colls);
        }
    } else {
        thrust::copy(host_pars.begin(), host_pars.end(), device_pars.begin());

        total_colls = 0;
        cudaMemcpyToSymbol(cnt, &total_colls, sizeof(total_colls));

        simulate_perf<<<num_blocks, num_threads>>>(host_n, host_l, device_pars_ptr, K, host_s + 0.0);

        cudaDeviceSynchronize();

        perf_resolve(host_pars_ptr, device_pars_ptr, host_s + 0.0);
    }

    print_statistics(host_s, host_pars_ptr);

    // cout << "collision count: " << totalP << " " << totalW << "\n";

    time(&end);

    double time_taken = end - start;
    // printf("Time lapsed: %10.8f\n", time_taken);
    
    device_pars_ptr = NULL;
    host_pars_ptr = NULL;
    cudaFree(colls);
    colls = NULL;
    device_pars.clear();
    device_pars.shrink_to_fit();
    host_pars.clear();
    host_pars.shrink_to_fit();
    

    return 0;
}
