#include <omp.h>
#include <bits/stdc++.h>

using namespace std;

const double EPS = 1e-8;
const double PI = acos(-1);
const double speedScale = 1e-3;

struct Point{
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

struct Particle {
    struct Point center = Point(0.0, 0.0);
    double radius;
    Velocity v = Velocity(0.0, 0.0);
    int id;

    int wColl;
    int pColl;

    double lastCol;

    Particle(struct Point C, double R, Velocity V, int ID) {
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
        if (newX >= limit - radius + EPS) newX = limit-radius;
        else if (newX + EPS <= radius) newX = radius;

        double newY = center.y + v.y * mtime;
        if (newY >= limit - radius + EPS) newY = limit-radius;
        else if (newY + EPS <= radius) newY = radius;

        center = Point(newX, newY);
    } //only check for wall "touching"

    bool Overlap(const Particle &oth) {
        Point shift = center - oth.center;
        double sqDist = shift * shift;
        double sqSumR = (radius + oth.radius) * (radius + oth.radius);
        return sqDist <= sqSumR - EPS; 
    }
};

struct Collision {
    vector<int> par;
    double ctime;
    int checkSum;

    Collision(vector<int> pars, double t) {
        par = pars;
        ctime = t;
    }

    bool operator< (const Collision &oth) const {
        if (ctime > oth.ctime + EPS)
            return false;
        else if (ctime + EPS < oth.ctime)
            return true;
        else {
            long l1 = par.size();
            long l2 = oth.par.size();
            if (l1 == l2 && l1 == 2)
                return par[0] < oth.par[0] || par[1] < oth.par[1];
            return par[0] < oth.par[0] || l1 < l2;
        }
    }
};

double CollideTimePar(Particle &a, Particle &b) {
    Velocity diffV = a.v - b.v;
    struct Point diffC = a.center - b.center;
    double distSq = diffC * diffC;
    double touchD = 4*a.radius*a.radius;
    if (distSq + EPS < touchD) {
		return -2.0;
        // fprintf(stderr, "Overlapping particles %d %d at %10.8f!!!\n", a.id, b.id, a.lastCol);
        // fprintf(stderr, "%d %10.8f %10.8f %10.8f %10.8f\n",a.id,a.center.x,a.center.y,a.v.x,a.v.y);
        // fprintf(stderr, "%d %10.8f %10.8f %10.8f %10.8f\n",b.id,b.center.x,b.center.y,b.v.x,b.v.y);
		// exit(1);
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
    /* return delta <= EPS 
            ? -2.0
            : (-bVal - sqrt(delta))/aVal >= EPS
                ? (-bVal - sqrt(delta))/aVal
                : (-bVal + sqrt(delta))/aVal;
    */
}

double CollideTimeWall(Particle &a, double boxLength) {
    double r = a.radius;
    double dist = boxLength - r;

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


struct Simulator{
    vector<Particle> par;
    double boxLength;

    bool AddParticle(const Particle &part, int scale){
        for(auto i : par) if (i.Overlap(part)) return 0; 

        ///Speed scale
        Particle newPar = part;
        newPar.v = newPar.v * (scale == -1 ? speedScale : scale);
        par.push_back(newPar);
        return 1;
    }

    bool isValidCollision (double ctime, double curTime) {
        return ctime >= curTime + EPS;
    }

    bool invalidWallCol (Collision &col, int ids) {
        return par[ids].pColl + par[ids].wColl != col.checkSum || par[ids].lastCol > col.ctime + EPS;
    }

    bool invalidParCol (Collision &col, int a, int b) {
        return par[a].pColl + par[a].wColl + par[b].pColl + par[b].wColl != col.checkSum 
                    || par[a].lastCol > col.ctime + EPS || par[b].lastCol > col.ctime + EPS;
    }

    bool wallColX (Collision &col, int ids) {
        return abs(par[ids].center.x - par[ids].radius) <= EPS || abs(boxLength - par[ids].center.x - par[ids].radius) <= EPS;
    }

    bool wallColY (Collision &col, int ids) {
        return abs(par[ids].center.y - par[ids].radius) <= EPS || abs(boxLength - par[ids].center.y - par[ids].radius) <= EPS;
    }

    bool wallNewValidTime(double t, int step) {
        return t + EPS > 0 && t <= step + EPS;
    }

    bool colNewValidTime (double t, int step) {
        return t > EPS && t <= step + EPS;
    }

    void parCol(int a, int b) {
        struct Point un = par[a].center - par[b].center;
        double unLen = sqrt(un * un);
        un = Point (un.x / unLen, un.y / unLen);
        struct Point ut = Point(-un.y, un.x);

        double v1n = un * par[a].v;
        double v1t = ut * par[a].v;
        double v2n = un * par[b].v;
        double v2t = ut * par[b].v;

        // db v1nf = v2n; db v2nf = v1n;
        par[a].v = Point(v2n * un.x + v1t * ut.x, v2n * un.y + v1t * ut.y);

        par[b].v = Point(v1n * un.x + v2t * ut.x, v1n * un.y + v2t * ut.y);

        par[a].pColl++;
        par[b].pColl++;
    }

    void wallCol(Collision &col, int p) {
        if (wallColX(col, p)) 
            par[p].v = Velocity(-par[p].v.x, par[p].v.y);
            
        if (wallColY(col, p)) 
            par[p].v = Velocity(par[p].v.x, -par[p].v.y);

        par[p].wColl++;
    }

    void StepSimulate(int step) {
        long n = par.size();
        for(long i = 0; i < n; i++) {
            Particle A = par[i];
            printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", 0, A.id, A.center.x, A.center.y, A.v.x, A.v.y);
        }

        for(int cnt = 1; cnt <= step; cnt++) {
            vector<Collision> event;
            int p0, p1;
            #pragma omp parallel for shared(event) private (p0,p1)
            for(p0 = 0; p0 < n; p0++) {
                for(p1 = 0; p1 < n; p1++) {
                    if (p1 == p0) continue;
                    
                    double colTime = CollideTimePar(par[p0], par[p1]);
                    if (colTime <= EPS || colTime > 1.0 + EPS) {
                        continue;
                    }
                    
                    vector<int> temp;
                    if (p0 < p1) {
                        temp.push_back(par[p0].id);
                        temp.push_back(par[p1].id);    
                    } else {
                        temp.push_back(par[p1].id);
                        temp.push_back(par[p0].id);
                    }
                   	
					#pragma omp critical 
                    event.push_back(Collision( temp, colTime ));
                }
                double wallTime = CollideTimeWall(par[p0], boxLength);
                if (wallTime + EPS < 0 || wallTime > 1.0 + EPS)
                    continue;
                
                vector<int> temp;
                temp.push_back(par[p0].id);
				#pragma omp critical
                event.push_back(Collision( temp, wallTime ));
            }
			#pragma omp barrier

            sort(event.begin(), event.end());

            vector<int> hasCollision (n, 0);
            for(Collision i : event) {
                if (i.par.size() == 1 && !hasCollision[i.par[0]]) {
                    ///wall collision
                    int ID = par[i.par[0]].id;
                    par[ID].Move(i.ctime, boxLength);
                    wallCol(i, ID);

                    hasCollision[ID] = 1;

                    if (abs(1.0 - i.ctime) > EPS)
                        par[ID].Move(1.0-i.ctime, boxLength);
                }
                else if (!hasCollision[i.par[0]] && !hasCollision[i.par[1]]) {
                    ///2 particles
 
                    int a = i.par[0];
                    int b = i.par[1];
                    par[a].Move(i.ctime, boxLength);
                    par[b].Move(i.ctime, boxLength);
                    parCol(a, b);

                    hasCollision[a] = 1;
                    hasCollision[b] = 1;

                    if (1.0 - i.ctime > EPS) {
                        par[a].Move(1.0-i.ctime, boxLength);
                        par[b].Move(1.0-i.ctime, boxLength);
                    }
                    
                }
            }
            
            int i;
            #pragma omp parallel for private(i)
            for(i = 0; i < n; i++) if (!hasCollision[i]) {
                par[i].Move(1.0, boxLength);
            }
			#pragma omp barrier
            
            if (cnt < step)
                for(i = 0; i < n; i++) {
                    Particle A = par[i];
                    printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", cnt, A.id, A.center.x, A.center.y, A.v.x, A.v.y);
                }
            else
                for(i = 0; i < n; i++) {
                    Particle A = par[i];
                    printf("%d %d %10.8f %10.8f %10.8f %10.8f %d %d\n", step,
                        A.id, A.center.x, A.center.y, A.v.x, A.v.y, A.pColl, A.wColl);
                }
        }
    }

    void addParCol(double colTime, int a, int b, priority_queue<Collision> *ev) {
        vector<int> temp;
        if (a < b) {
            temp.push_back(par[a].id);
            temp.push_back(par[b].id);    
        } else {
            temp.push_back(par[b].id);
            temp.push_back(par[a].id);
        }
        Collision col (temp, -colTime);
        col.checkSum = par[a].wColl + par[a].pColl + par[b].wColl + par[b].pColl;
        (*ev).push(col);
        
    }

    void addWallCol(double wallTime, int p, priority_queue<Collision> *ev) {
        vector<int> temp = {par[p].id};
        Collision col (temp, -wallTime);
        col.checkSum = par[p].wColl + par[p].pColl;
        (*ev).push(col);
    }

    void FastSimulate(int step) {
        ///TODO
        
        priority_queue<Collision> ev;
        int p0, p1;
        int n = par.size();
		for(int i = 0; i < n; i++) {
            Particle A = par[i];
            printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", 0, A.id, A.center.x, A.center.y, A.v.x, A.v.y);
        }
        #pragma omp parallel for shared(ev) private (p0,p1)
        for(p0 = 0; p0 < n; p0++) {
            for(p1 = 0; p1 < n; p1++) {
                if (p1 == p0) continue;
                
                double colTime = CollideTimePar(par[p0], par[p1]);

                if (!colNewValidTime(colTime, step)) 
                    continue;
                
				#pragma omp critical
                addParCol(colTime, p0, p1, &ev);
            }
            double wallTime = CollideTimeWall(par[p0], boxLength);

            if (!wallNewValidTime(wallTime, step))
                continue;
                  
			#pragma omp critical
            addWallCol(wallTime, p0, &ev);
        }
		#pragma omp barrier

        double curTime = 0.0;
        double lapse;
        while (!ev.empty()) {
            Collision col = ev.top();
			col.ctime = -col.ctime;
            ev.pop();
            if (!isValidCollision(col.ctime, curTime)) continue;
            lapse = col.ctime - curTime;
            curTime = col.ctime;

			int i;
            #pragma omp parallel for private (i)
            for (i = 0; i < n; i++) 
                par[i].Move(lapse, boxLength);
			#pragma omp barrier

            if (col.par.size() == 1) {
                // wall
                int ids = col.par[0];
                if (invalidWallCol(col, ids))
                    continue;
                
                wallCol(col, ids);
                par[ids].lastCol = ceil(col.ctime);

                double wallTime = CollideTimeWall(par[ids], boxLength) + curTime;
                
                if (wallNewValidTime(wallTime, step)) {
                    if (wallTime <= par[ids].lastCol + EPS)
						wallTime = par[ids].lastCol;
					addWallCol(wallTime, ids, &ev);
                } 
                
                #pragma omp parallel for private (i) shared (ev)
                for (i = 0; i < n; i++) {
                    if (i == ids) continue;
                    double colTime = CollideTimePar(par[i], par[ids]) + curTime;

                    if (!colNewValidTime(colTime, step) || colTime <= par[ids].lastCol + EPS) {
                        continue;
                    }

					#pragma omp critical
                    addParCol(colTime, i, ids, &ev);
                }
            } else if (col.par.size() == 2){
                // 2 pars
                int a = col.par[0];
                int b = col.par[1];
                if (invalidParCol(col, a, b)) 
                    continue;

                parCol(a, b);

                par[a].lastCol = ceil(col.ctime);
                par[b].lastCol = ceil(col.ctime);
				
				
				double wallA = CollideTimeWall(par[a], boxLength) + curTime;
                
                if (wallNewValidTime(wallA, step)) {
                    if (wallA <= par[a].lastCol + EPS)
                        wallA = par[a].lastCol;
                    addWallCol(wallA, a, &ev);
                }
				
				double wallB = CollideTimeWall(par[b], boxLength) + curTime;

                if (wallNewValidTime(wallB, step)) {
                    if (wallB <= par[b].lastCol + EPS)
                        wallB = par[b].lastCol;
                    addWallCol(wallB, b, &ev);
                }
                int i;
                
				#pragma omp parallel for shared(ev) private (i)                
                for (i = 0; i < n; i++) {
                    if (i != a) {
                        double colTime = CollideTimePar(par[i], par[a]) + curTime;

                        if (colNewValidTime(colTime, step) && colTime > par[a].lastCol + EPS) {    
                            #pragma omp critical
							addParCol(colTime, i, a, &ev);
                        }
                    }
                    if (i != b) {
                        double colTime = CollideTimePar(par[i], par[b]) + curTime;

                        if (colNewValidTime(colTime, step) && colTime > par[b].lastCol + EPS) {
                            #pragma omp critical
							addParCol(colTime, i, b, &ev);
                        }
                    }
            	}
            }
        }

        if (step - curTime > EPS) {
            int k;
            #pragma omp parallel for private(k)
            for (k = 0; k < n; k++)
                par[k].Move(step - curTime, boxLength);
        	#pragma omp barrier
		}

        for (int i = 0; i < n; i++) {
            Particle A = par[i];
            printf("%d %d %10.8f %10.8f %10.8f %10.8f %d %d\n", step,
                A.id, A.center.x, A.center.y, A.v.x, A.v.y, A.pColl, A.wColl);
        }
    }
} mySim;

double rng01(mt19937 &myRng) {
    unsigned value = myRng();
    return (value + .0) / myRng.max();
}

Particle RandomParticle(double boxLength, double radius, int index, 
    mt19937& myRng) {
    
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

int main(int argc, char *argv[]) {
    int N, S;
    double L, r;
    string progType;
    cin >> N >> L >> r >> S >> progType; 

    mySim.boxLength = L;

    long idx;
    double x, y, vx, vy;
    if (scanf("%ld %lf %lf %lf %lf", &idx, &x, &y, &vx, &vy) != EOF) {
        // cout << "read from input\n";
        long i = 0;
        mySim.AddParticle(Particle(Point(x,y),
                    r,
                    Velocity(vx, vy),
                    idx), 1);
        for(i = 1; i < N; i++) {
            
            cin >> idx >> x >> y >> vx >> vy;
            mySim.AddParticle(Particle(Point(x,y),
                    r,
                    Velocity(vx, vy),
                    idx), 1);
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
                if (!mySim.AddParticle(newP, -1)) i--;
            }
        } else {
            mt19937 myRng(time(NULL));

            for(int i = 0; i < N; i++) {
                Particle newP = RandomParticle(L, r, i, myRng); 
                if (!mySim.AddParticle(newP, -1)) i--;
            }
        }   
    }
    

    /*
    Change the number to run in a different level of parallelism
    */ 
    omp_set_num_threads(16);


    if (progType == "print") {
        mySim.StepSimulate(S);
    } else {
        mySim.FastSimulate(S);
    }

    return 0;
}
