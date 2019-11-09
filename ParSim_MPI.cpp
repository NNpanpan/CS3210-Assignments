#include <stdio.h>
#include <string.h>
#include <bits/stdc++.h>

using namespace std;

const double eps = 1e-8;
const double PI = acos(-1);

long N, S;
double L, r;
/*
-----TRUCTURES-----
*/
struct Point
{ 
    double x, y;

    Point(double cx, double cy) {
        x = cx;
        y = cy;
    }

    Point operator + (const Point &oth) {
        return Point(x + oth.x, y + oth.y);
    }
    Point operator - (const Point &oth) {
        return Point(x - oth.x, y - oth.y);
    }
    Point operator * (double k) {
        return Point(x*k, y*k);
    }
    double operator * (const Point &oth) {
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

    Particle(struct Point C, double R, Velocity V, long ID) {
        center = C;
        radius = R;
        v = V;
        id = ID;
        wColl = 0;
        pColl = 0;
        lastCol = -1.0;
    }

    void Move(double mtime, double limit){
        double newX = center.x + v.x * mtime;
        if (newX >= limit - radius + eps) newX = limit-radius;
        else if (newX + eps <= radius) newX = radius;
        
        double newY = center.y + v.y * mtime;
        if (newY >= limit - radius + eps) newY = limit-radius;
        else if (newY + eps <= radius) newY = radius;

        center = Point(newX, newY);
    } //only check for wall "touching"

    bool Overlap(const Particle &oth) {
        Point shift = center - oth.center;
        double sqDist = shift * shift;
        double squmR = (radius + oth.radius) * (radius + oth.radius);
        return sqDist <= squmR - eps; 
    }
};

struct Collision 
{
    long par1, par2;
    double ctime;
    long checkum;

    Collision(long p1, long p2, double t, long sum = 0) {
        par1 = p1 > p2 ? p1 : p2;
        par2 = p1 + p2 - par1;
        ctime = t;
        checkum = sum;
    }

    bool operator < (const Collision &oth) const {
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

Collision* colls;
vector<Particle> pars;

struct collHeap 
{
    long n;

    collHeap() {
        n = 0;
    }

    void shiftUp(long index) {
        long k = index;
        while (k > 1 && colls[k] < colls[k/2]) {
            Collision temp = colls[k];
            colls[k] = colls[k/2];
            colls[k/2] = temp;
            k /= 2;
        }
    }

    void shiftDown(long index) {
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

    void insert(Collision C) {
        n++;
        colls[n] = C;
        shiftUp(n);
    }

    Collision pop() {
        Collision ret = colls[1];
        colls[1] = colls[n];
        n--;
        shiftDown(1);
        return ret;
    }

    bool isEmpty() {
        return n == 0;
    }

    void build() {
        for (long k = n/2; k >= 1; k--) {
            shiftDown(k);
        }
    }
};


/*
-----COLLIION TIME-----
*/
double CollideTimePar(Particle &a, Particle &b) {
	Velocity diffV = a.v - b.v;
    struct Point diffC = a.center - b.center;
    double distq = diffC * diffC;
    double touchD = 4*a.radius*a.radius;
    if (distq + eps < touchD) {
		return 0.0;
    }

    double bVal = diffV * diffC;
    double aVal = diffV * diffV;
    double delta = bVal*bVal - aVal*(distq - touchD);
    if (delta <= eps) {
        return -2.0;
    } else {
        double ret1 = (-bVal - sqrt(delta))/aVal;
        double ret2 = (-bVal + sqrt(delta))/aVal;
        return ret1 >= 0.0 + eps ? ret1 : ret2;
    }
}

double CollideTimeParPerf(Particle &a, Particle &b, double cur) {
	Velocity diffV = a.v - b.v;
    struct Point diffC = a.center - b.center;
    double distq = diffC * diffC;
    double touchD = 4*a.radius*a.radius;
    if (distq + eps < touchD) {
		return -2.0;
    }

    double bVal = diffV * diffC;
    double aVal = diffV * diffV;
    double delta = bVal*bVal - aVal*(distq - touchD);
    if (delta <= eps) {
        return -2.0;
    } else {
        double ret1 = (-bVal - sqrt(delta))/aVal + cur;
        double ret2 = (-bVal + sqrt(delta))/aVal + cur;
        double maxLastCol = max (a.lastCol, b.lastCol);
        if (ret1 + eps < maxLastCol && maxLastCol + eps < ret2) return maxLastCol;
        else if (maxLastCol + eps < ret1) return ret1;
        else return -2.0;
    }
}

double CollideTimeWall(Particle &a, double boxLength) {
    double r = a.radius;

    double xTime, yTime;
    if (a.v.x > 0.0 + eps) {
        xTime = abs(boxLength - r - a.center.x)/a.v.x;
    } else if (a.v.x <= -eps) {
        xTime = abs((a.center.x - r)/a.v.x);
    } else {
        xTime = DBL_MAX;
    }

    if (a.v.y > 0.0 + eps) {
        yTime = abs(boxLength - r - a.center.y)/a.v.y;
    } else if (a.v.y <= -eps) {
        yTime = abs((a.center.y - r)/a.v.y);
    } else {
        yTime = DBL_MAX;
    }

    return xTime > yTime + eps ? yTime : xTime;
}



/*
-----CONDITION TESTING-----
*/
bool isValidCollisionTime (double ctime, double curTime, double lim) 
{
    return ctime + eps > curTime && ctime + eps < lim;
}

bool isValidChecksum (Particle* pars, long id1, long id2, long checkum) 
{
    if (id2 == -1) {
        return pars[id1].wColl + pars[id1].pColl == checkum;
    } else {
        return pars[id1].wColl + pars[id1].pColl
            + pars[id2].wColl + pars[id2].pColl == checkum;
    }
}

bool validLastCols (Particle* pars, long id1, long id2, double ctime) 
{
    if (id2 == -1) {
        return pars[id1].lastCol <= ctime + eps;
    } else {
        return pars[id1].lastCol <= ctime + eps && pars[id2].lastCol <= ctime + eps;
    }
}

bool wallColX (int ids, Particle* pars) 
{
    return abs(pars[ids].center.x - pars[ids].radius) <= eps 
        || abs(L - pars[ids].center.x - pars[ids].radius) <= eps;
}

bool wallColY (int ids, Particle* pars) 
{
    return abs(pars[ids].center.y - pars[ids].radius) <= eps 
        || abs(L - pars[ids].center.y - pars[ids].radius) <= eps;
}



/*
-----COLLIION RESOLVING-----
*/
void parCol(Particle* pars, int a, int b) 
{
    struct Point un = pars[a].center - pars[b].center;
    double unLen = sqrt(un * un);
    un = Point (un.x / unLen, un.y / unLen);
    struct Point ut = Point(-un.y, un.x);

    double v1n = un * pars[a].v;
    double v1t = ut * pars[a].v;
    double v2n = un * pars[b].v;
    double v2t = ut * pars[b].v;

    // db v1nf = v2n; db v2nf = v1n;
    pars[a].v = Point(v2n * un.x + v1t * ut.x, v2n * un.y + v1t * ut.y);

    pars[b].v = Point(v1n * un.x + v2t * ut.x, v1n * un.y + v2t * ut.y);

    pars[a].pColl++;
    pars[b].pColl++;
}

void wallCol(Particle* pars, int p) 
{
    if (wallColX(p, pars)) 
        pars[p].v = Velocity(-pars[p].v.x, pars[p].v.y);
        
    if (wallColY(p, pars)) 
        pars[p].v = Velocity(pars[p].v.x, -pars[p].v.y);

    pars[p].wColl++;
}

long totalP = 0; // to count total collisions between particles
long totalW = 0; // to count total wall collisions
/*
-----PRINT MODE-----
*/
void simulate_step(long n, double l, Particle* device_pars, long lim)
{
	long len = n;
    long i = -1; 	
    long k = -1;    
	//printf("kernel on thread: %ld\n", i+k);
    for (long cnt1 = 0; i + cnt1 < len && cnt1 < lim; cnt1 ++) {
        for (long cnt2 = 0; k + cnt2 < len && cnt2 < lim; cnt2 ++) {
			// if (i+cnt1 < k+cnt2) printf("particles %ld and %ld\n", i+cnt1, k+cnt2);
            if (i + cnt1 == k + cnt2) {
                double wallTime = CollideTimeWall(device_pars[i + cnt1], l);
                
				if (wallTime + eps < 0 || wallTime > 1.0 + eps)
                    continue;
                // printf("%ld will hit wall at %10.8f\n", i+cnt1, wallTime);
                // add to heap
                // printf("collision %ld\n", j);
            } else {
                long p1 = i + cnt1;
                long p2 = k + cnt2;
                double colTime = CollideTimePar(device_pars[p1], device_pars[p2]);
                if (colTime + eps < 0 || colTime > 1.0 + eps) {
                    continue;
                }
                
                
                // printf("collision %ld\n", j);
                // add to heap
            }
        }
    }
    
}

void step_resolve(Particle* pars, Particle* device_pars, long coll_size) 
{
    // cout << "sorting size " << coll_size << "\n";
    std::sort(colls, colls + coll_size);
    vector<int> hasCollision (N, 0);
	
	// cout << "resolving\n";
    for(long i = 0; i < coll_size; i++) {
        Collision c = colls[i];
        if (c.par2 == -1 && !hasCollision[c.par1]) {
            ///wall collision

			long ID = c.par1;
			// printf("particle %ld hit wall at %10.8f\n", ID, c.ctime);
            pars[ID].Move(c.ctime, L);
            wallCol(pars, ID);

            hasCollision[ID] = 1;

            if (1.0 - c.ctime > eps)
                pars[ID].Move(1.0-c.ctime, L);

            totalW++;
        }
        else if (!hasCollision[c.par1] && !hasCollision[c.par2]) {
            ///2 particles

            long a = c.par1;
            long b = c.par2;
			// printf("particles %ld and %ld hit at %10.8f\n", a, b, c.ctime);
            pars[a].Move(c.ctime, L);
            pars[b].Move(c.ctime, L);
            parCol(pars, a, b);

            hasCollision[a] = 1;
            hasCollision[b] = 1;

            if (1.0 - c.ctime > eps) {
                pars[a].Move(1.0-c.ctime, L);
                pars[b].Move(1.0-c.ctime, L);
            }

            totalP++;
        }
    }
    // cout << "done resolving\n";
    for(long i = 0; i < N; i++) if (!hasCollision[i]) {
        pars[i].Move(1.0, L);
    }
    // cout << "done moving. new loop ready\n";
}

/*
-----PERF MODE-----
*/
void simulate_perf(long n, long l, Particle* device_pars, long lim, double s) 
{
    long len = n;
    long i = lim;
    long k = lim;
      
    for (long cnt1 = 0; i + cnt1 < len && cnt1 < lim; cnt1 ++) {
        for (long cnt2 = 0; k + cnt2 < len && cnt2 < lim; cnt2 ++) {
            // if (i+cnt1 < k+cnt2) printf("particles %ld and %ld\n", i+cnt1, k+cnt2);
            if (i + cnt1 == k + cnt2) {
                long p = i+cnt1;
                double wallTime = CollideTimeWall(device_pars[p], l);
                
				if (wallTime + eps < 0 || wallTime > s + eps)
                    continue;
                
                // printf("collision %ld\n", pos);
                // add to heap
            } else {
                long p1 = i + cnt1;
                long p2 = k + cnt2;
                if (p1 > p2) continue;
                // printf("particles %ld and %ld\n", p1, p2);

                double colTime = CollideTimeParPerf(device_pars[p1], device_pars[p2], 0.0);
                if (colTime + eps < 0 || colTime > s + eps) {
                    continue;
                }
                
                // printf("collision %ld\n", pos);
                // add to heap
            }
        }
    }
}

 void parallelRecalcingle(Particle* pars, long id, double step, long n, long l, long K, double cur)
{
    long i =  K;
    for (long ct = 0; i + ct < n && ct < K; ct++) {
        if (i + ct == id) {
            double wallTime = CollideTimeWall(pars[id], l) + cur;
                
            if (wallTime + eps < cur || wallTime > step + eps)
                continue;
            
            if (wallTime <= pars[id].lastCol + eps)
                wallTime = pars[id].lastCol;

            long sum = pars[id].wColl + pars[id].pColl;
            // add to heap
        } else {
            double colTime = CollideTimeParPerf(pars[i+ct], pars[id], cur);
            if (colTime + eps < cur || colTime > step + eps) {
                continue;
            }
            
            long sum = pars[id].wColl + pars[id].pColl
                    + pars[i+ct].wColl + pars[i+ct].pColl;
            // add to heap
        }
    }
}

 void parallelRecalcDouble(Particle* pars, long id1, long id2, double step, long n, long l, long K, double cur)
{
    long i = K;
    for (long ct = 0; i + ct < n && ct < K; ct++) {
        long p = i + ct;
        if (p == id1) {
            double wallTime = CollideTimeWall(pars[id1], l) + cur;
                
            if (wallTime + eps < cur || wallTime > step + eps)
                continue;
            
            if (wallTime <= pars[id1].lastCol + eps)
                wallTime = pars[id1].lastCol;

            long sum = pars[id1].wColl + pars[id1].pColl;
        } else if (p == id2) {
            double wallTime = CollideTimeWall(pars[id2], l) + cur;
                
            if (wallTime + eps < cur || wallTime > step + eps)
                continue;

            if (wallTime <= pars[id2].lastCol + eps)
                wallTime = pars[id2].lastCol;
            
            long sum = pars[id2].wColl + pars[id2].pColl;
        } else {    
            double colTime1 = CollideTimeParPerf(pars[p], pars[id1], cur);
            double colTime2 = CollideTimeParPerf(pars[p], pars[id2], cur);
            
            if (colTime1 + eps >= cur && colTime1 <= step + eps) {
                long sum = pars[id1].wColl + pars[id1].pColl
                    + pars[p].wColl + pars[p].pColl;
            }

            if (colTime2 + eps >= cur && colTime2 <= step + eps) {
                long sum = pars[id2].wColl + pars[id2].pColl
                    + pars[p].wColl + pars[p].pColl;
            }
            
        }
    }
}

 void perf_resolve(Particle* pars, Particle* device_pars_ptr, double step) 
{
}



/*
-----TATISTICS-----
*/
 void print_particles(int step)
{
    int i;
    for (i = 0; i < N; i++) {
        Particle A = pars[i];
        printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", step, i, A.center.x, A.center.y,
            A.v.x, A.v.y);
    }
}

 void print_statistics(int num_step)
{
    int i;
    for (i = 0; i < N; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f %ld %ld\n", num_step, i, pars[i].center.x,
            pars[i].center.y, pars[i].v.x, pars[i].v.y,
            pars[i].pColl, pars[i].wColl);
    }
}


/*
-----RANDOMIZER-----
*/
 double rng01(mt19937 &myRng) 
{
    unsigned value = myRng();
    return (value + .0) / myRng.max();
}

 Particle RandomParticle(double boxLength, double radius, int index, 
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

 bool AddParticle(const Particle &part, int pos)
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

    cin >> N >> L >> r >> S >> mode;

    long idx;
    double x, y, vx, vy;
    if (scanf("%ld %lf %lf %lf %lf", &idx, &x, &y, &vx, &vy) != EOF) {
        // cout << "read from input\n";
        long i = 0;
        pars[i] = (Particle(Point(x,y),
                    r,
                    Velocity(vx, vy),
                    idx));
        for(i = 1; i < N; i++) {
            
            cin >> idx >> x >> y >> vx >> vy;
            pars[i] = (Particle(Point(x,y),
                    r,
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
            for(int i = 0; i < N; i++) {
                Particle newP = RandomParticle(L, r, i, myRng); 
                if (!AddParticle(newP, i)) i--;
            }
        } else {
            mt19937 myRng(time(NULL));

            for(int i = 0; i < N; i++) {
                Particle newP = RandomParticle(L, r, i, myRng); 
                if (!AddParticle(newP, i)) i--;
            }
        }   
    }



    long total_colls;
    
    print_particles(0);

    if (mode == "print") {
        // cout << "print mode\n";
        
        for (long step = 0; step < S; step++) {
			// cout << "step: " << step << "\n";
            if (step > 0) print_particles(step);

        }
    } else {
    }

    print_statistics(S);

    // cout << totalP << " " << totalW << "\n";
    

    return 0;
}
