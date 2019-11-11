#include <stdio.h>
#include <string.h>
#include <bits/stdc++.h>
#include <mpi.h>

// using namespace std;

const double eps = 1e-8;
const double PI = acos(-1);

int N, S;
double L, r;
std::string mode;

int rank, size;

/*
-----STRUCTURES-----
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
    int id;

    int wColl;
    int pColl;


    Particle(struct Point C, double R, Velocity V, int ID) {
        center = C;
        radius = R;
        v = V;
        id = ID;
        wColl = 0;
        pColl = 0;
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
    int par1, par2;
    double ctime;

    Collision(int p1, long p2, double t) {
        par1 = p1 > p2 ? p2 : p1;
        par2 = p1 + p2 - par1;
        ctime = t;
    }

    bool operator < (const Collision &oth) const {
        if (ctime > oth.ctime + eps)
            return false;
        else if (ctime + eps < oth.ctime)
            return true;
        else {
            if (par1 == oth.par1) return par2 < oth.par2;
            return par1 < oth.par1;
        }
    }
};

// Necessary types
MPI_Datatype mpi_point_type;
MPI_Datatype mpi_velocity_type;
MPI_Datatype mpi_particle_type;
MPI_Datatype mpi_collision_type;

// Containers
Collision* colls;
Particle* pars;

struct collHeap 
{
    int n;

    collHeap() {
        n = 0;
    }

    void shiftUp(int index) {
        int k = index;
        while (k > 1 && colls[k] < colls[k/2]) {
            Collision temp = colls[k];
            colls[k] = colls[k/2];
            colls[k/2] = temp;
            k /= 2;
        }
    }

    void shiftDown(int index) {
        int k = index;
        while (k*2 <= n) {
            int dir = k * 2;
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

    void buid() {
        for (int k = n/2; k >= 1; k--) {
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

int totalP = 0; // to count total collisions between particles
int totalW = 0; // to count total wall collisions

/*
-----RESOLVING COLLISIONS-----
*/
void step_resolve(int coll_size) 
{
    // cout << "sorting size " << coll_size << "\n";
    std::sort(colls, colls + coll_size);
    std::vector<int> hasCollision (N, 0);
	
	// cout << "resolving\n";
    for(int i = 0; i < coll_size; i++) {
        Collision c = colls[i];
        if (c.par2 == -1 && !hasCollision[c.par1]) {
            ///wall collision

			int ID = c.par1;
			// printf("particle %d hit wall at %10.8f\n", ID, c.ctime);
            pars[ID].Move(c.ctime, L);
            wallCol(pars, ID);

            hasCollision[ID] = 1;

            if (1.0 - c.ctime > eps)
                pars[ID].Move(1.0-c.ctime, L);

            totalW++;
        }
        else if (!hasCollision[c.par1] && !hasCollision[c.par2]) {
            ///2 particles

            int a = c.par1;
            int b = c.par2;
			// printf("particles %d and %ld hit at %10.8f\n", a, b, c.ctime);
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
    for(int i = 0; i < N; i++) if (!hasCollision[i]) {
        pars[i].Move(1.0, L);
    }
    // cout << "done moving. new loop ready\n";
}


/*
-----STATISTICS-----
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

 void print_statistics()
{
    int i;
    for (i = 0; i < N; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f %d %d\n", S, i, pars[i].center.x,
            pars[i].center.y, pars[i].v.x, pars[i].v.y,
            pars[i].pColl, pars[i].wColl);
    }
}


/*
-----RANDOMIZER-----
*/
 double rng01(std::mt19937 &myRng) 
{
    unsigned value = myRng();
    return (value + .0) / myRng.max();
}

 Particle RandomParticle(double boxLength, double radius, int index, 
    std::mt19937& myRng) 
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

/*
-----INPUT-----
*/
void readInput(int argc, char** argv) {
    std::cin >> N >> L >> r >> S >> mode;

    pars = (Particle*) malloc(sizeof(Particle) * N);

    int idx;
    double x, y, vx, vy;
    if (scanf("%d %lf %lf %lf %lf", &idx, &x, &y, &vx, &vy) != EOF) {
        // cout << "read from input\n";
        int i = 0;
        pars[i] = (Particle(Point(x,y),
                    r,
                    Velocity(vx, vy),
                    idx));
        for(i = 1; i < N; i++) {
            std::cin >> idx >> x >> y >> vx >> vy;
            pars[i] = (Particle(Point(x,y),
                    r,
                    Velocity(vx, vy),
                    idx));
        }
    } else {
		// cout << "randomize\n";
        std::string seed = "agony";
        if (argc == 2) {
            seed = std::string(argv[1]);
            std::seed_seq ss (seed.begin(), seed.end());
            std::mt19937 myRng(ss);
            for(int i = 0; i < N; i++) {
                Particle newP = RandomParticle(L, r, i, myRng); 
                if (!AddParticle(newP, i)) i--;
            }
        } else {
            std::mt19937 myRng(time(NULL));

            for(int i = 0; i < N; i++) {
                Particle newP = RandomParticle(L, r, i, myRng); 
                if (!AddParticle(newP, i)) i--;
            }
        }
    }
}

/*
-----SETUP MPI TYPES-----
*/
void type_setup() {
    const int pnt_items = 2;
    const int par_items = 6;
    const int col_items = 3;

    int block_length_pnt[2] = {1,1};
    int block_length_par[6] = {1,1,1,1,1,1};
    int block_length_col[3] = {1,1,1};

    //Setup for Point and Velocity
    MPI_Datatype types_pnt[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets_pnt[2];
    offsets_pnt[0] = offsetof(Point, x);
    offsets_pnt[1] = offsetof(Point, y);
    MPI_Type_create_struct(pnt_items, block_length_pnt, offsets_pnt, types_pnt, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);
    MPI_Type_create_struct(pnt_items, block_length_pnt, offsets_pnt, types_pnt, &mpi_velocity_type);
    MPI_Type_commit(&mpi_velocity_type);

    //Setup for Particle
    MPI_Datatype types_par[6] = {mpi_point_type, MPI_DOUBLE, mpi_velocity_type, MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint offsets_par[6];
    offsets_par[0] = offsetof(Particle, center);
    offsets_par[1] = offsetof(Particle, radius);
    offsets_par[2] = offsetof(Particle, v);
    offsets_par[3] = offsetof(Particle, id);
    offsets_par[4] = offsetof(Particle, wColl);
    offsets_par[5] = offsetof(Particle, pColl);
    MPI_Type_create_struct(par_items, block_length_par, offsets_par, types_par, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);

    //Setup for Collision
    MPI_Datatype types_col[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
    MPI_Aint offsets_col[3];
    offsets_col[0] = offsetof(Collision, par1);
    offsets_col[1] = offsetof(Collision, par2);
    offsets_col[2] = offsetof(Collision, ctime);
    MPI_Type_create_struct(col_items, block_length_col, offsets_col, types_col, &mpi_collision_type);
    MPI_Type_commit(&mpi_collision_type);
}

/*
-----MASTER AND WORKER-----
*/
void master() {
    int NS_values[2] = {N, S};
    double Lr_values[2] = {L, r};
    MPI_Bcast(NS_values, 2, MPI_INT, size-1, MPI_COMM_WORLD);
    MPI_Bcast(Lr_values, 2, MPI_DOUBLE, size-1, MPI_COMM_WORLD);
    
    MPI_Status status;
    colls = (Collision*) malloc(sizeof(Collision) * N*N);
    for (int step = 1; step <= S; step++) {
        MPI_Bcast(pars, N, mpi_particle_type, size-1, MPI_COMM_WORLD);
        int num_of_colls = 0; int cnt = 0;
        for (int wrkr = 0; wrkr < size-1; wrkr++) {
            MPI_Recv(&cnt, 1, MPI_INT, wrkr, NULL, MPI_COMM_WORLD, &status);
            MPI_Recv(colls + num_of_colls, cnt, mpi_collision_type, wrkr, NULL, MPI_COMM_WORLD, &status);
            num_of_colls += cnt;
        }

        step_resolve(num_of_colls);

        if (mode == "print" && step < S)
            print_particles(step);

        if (step == S)
            print_statistics();
    }
}

void worker() {
    int NS_values[2];
    double Lr_values[2];
    MPI_Bcast(NS_values, 2, MPI_INT, size-1, MPI_COMM_WORLD);
    MPI_Bcast(Lr_values, 2, MPI_DOUBLE, size-1, MPI_COMM_WORLD);
    N = NS_values[0];
    S = NS_values[1];
    L = Lr_values[0];
    r = Lr_values[1];

    pars = (Particle*) malloc(sizeof(Particle) * N);
    colls = (Collision*) malloc(sizeof(Collision) * (N*N+1)/size);
    for (int step = 1; step <= S; step++) {
        MPI_Bcast(pars, N, mpi_particle_type, size-1, MPI_COMM_WORLD);
        int cnt = 0;

        for (int i = 0; i < (N + size - 2) / (size - 1); i++) {
            int base = (i % 2) ? (i + 1) * (size - 1) - rank - 1 : i * (size - 1) + rank;
            if (base < N) {
                for (int k = base + 1; k < N; k++) {
                    double ctime = CollideTimePar(pars[base], pars[k]);
                    if (ctime + eps >= 0 && ctime <= 1.0 + eps) {
                        colls[cnt] = Collision(base, k, ctime);
                        cnt++;
                    }
                }
                double wTime = CollideTimeWall(pars[base], L);
                if (wTime + eps >= 0 && wTime <= 1.0 + eps) {
                    colls[cnt] = Collision(-1, base, wTime);
                    cnt++;
                }
            }
        }

        MPI_Send(&cnt, 1, MPI_INT, size-1, NULL, MPI_COMM_WORLD);
        MPI_Send(colls, cnt, MPI_INT, size-1, NULL, MPI_COMM_WORLD);
    }
    
}


int main(int argc, char** argv)
{

    MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    type_setup();

    if (rank == size - 1) {// master process function
        readInput(argc, argv);
        master();
    } else {// worker process function
        worker();
    }

    free(pars);
    free(colls);    

    return 0;
}
