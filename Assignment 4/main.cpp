#include <fstream>
#include <sstream>
#include "voronoi.h"
#include <map>
int debug = 0;
//Statements for debugging purposes
void update_k(){
    ++debug;
}
void print_k(){
    cerr << "K is :" << debug << endl;
}
int VoronoiPointCompare(const void *p1, const void *p2)
{
    VoronoiPoint *s1 = (VoronoiPoint*)p1, *s2 = (VoronoiPoint*)p2;
    update_k();
    if (s1->y < s2->y) return -1;
    if (s1->y > s2->y) return 1;
    update_k();
    if (s1->x < s2->x) return -1;
    if (s1->x > s2->x) return 1;
    update_k();
    return 0;
}

VoronoiPoint::VoronoiPoint(double nx, double ny)
{
    x = nx;
    y = ny;
    update_k();
}

VoronoiPoint::VoronoiPoint()
{
    x = 0.0;
    y = 0.0;
}

Voronoi::Voronoi()
{
    siteidx = 0;
    sites = 0;
    allMemoryList = new FreeNodeArrayList;
    allMemoryList->memory = 0;
    allMemoryList->next = 0;
    update_k();
    currentMemoryBlock = allMemoryList;
    allEdges = 0;
    update_k();
    iteratorEdges = 0;
    minDistanceBetweenSites = 0;
}

Voronoi::~Voronoi()
{
}

vector<VEdge> Voronoi::ComputeVoronoiGraph(vector<VoronoiPoint*> p, double minY, double maxY)
{
    iteratorEdges = allEdges;
    cleanup();
    cleanupEdges();
    update_k();
    int i;

    minDistanceBetweenSites = 0;
    update_k();
    nsites = p.size();
    plot = 0;
    update_k();
    triangulate = 0;
    debug = 1;
    update_k();
    sorted = 0;
    freeinit(&sfl, sizeof(Site));
    update_k();

    sites = (Site *)myalloc(nsites*sizeof(*sites));
    xmin = p[0]->x;
    update_k();
    ymin = p[0]->y;
    xmax = p[0]->x;
    update_k();
    ymax = p[0]->y;
    update_k();
    for (i = 0; i< nsites; i++)
    {
        sites[i].coord.x = p[i]->x;
        update_k();
        sites[i].coord.y = p[i]->y;
        update_k();
        sites[i].sitenbr = i;
        update_k();
        sites[i].refcnt = 0;
        update_k();
        if (p[i]->x < xmin)
            xmin = p[i]->x;
        else if (p[i]->x > xmax)
            xmax = p[i]->x;
        update_k();
        if (p[i]->y < ymin)
            ymin = p[i]->y;
        else if (p[i]->y > ymax)
            ymax = p[i]->y;
    }
    update_k();
    qsort(sites, nsites, sizeof(*sites), VoronoiPointCompare);
    update_k();
    siteidx = 0;
    geominit();
    double temp = 0;
    update_k();
    if (minY > maxY)
    {
        temp = minY;
        minY = maxY;
        update_k();
        maxY = temp;
    }
    if (minY > maxY)
    {
        temp = minY;
        update_k();
        minY = maxY;
        maxY = temp;
        update_k();
    }
    update_k();
    borderMinX = minY;
    update_k();
    borderMinY = minY;
    borderMaxX = maxY;
    update_k();
    borderMaxY = maxY;
    update_k();
    siteidx = 0;
    voronoi(triangulate);
    update_k();
    p.clear();
    cleanup();
    cleanupEdges();
    clean();
    return total_edges;

}
void Voronoi::clean()
{
    delete[] sites;
        update_k();
    delete[] PQhash;
        update_k();
    delete[] currentMemoryBlock;
        update_k();

    delete[] allEdges;
        update_k();

    delete[] iteratorEdges;
        update_k();

    delete[] ELhash;
}

bool Voronoi::ELinitialize()
{
    int i;
    freeinit(&hfl, sizeof **ELhash);
        update_k();

    ELhashsize = 2 * sqrt_nsites;
        update_k();

    ELhash = (Halfedge **)myalloc(sizeof *ELhash * ELhashsize);
            update_k();

    if (ELhash == 0)
        return false;
              update_k();

    for (i = 0; i<ELhashsize; i += 1) ELhash[i] = (Halfedge *)NULL;
    ELleftend = HEcreate((Edge *)NULL, 0);
    update_k();

    ELrightend = HEcreate((Edge *)NULL, 0);
        update_k();

    ELleftend->ELleft = (Halfedge *)NULL;
        update_k();

    ELleftend->ELright = ELrightend;
        update_k();

    ELrightend->ELleft = ELleftend;
        update_k();

    ELrightend->ELright = (Halfedge *)NULL;
        update_k();

    ELhash[0] = ELleftend;
        update_k();

    ELhash[ELhashsize - 1] = ELrightend;
        update_k();


    return true;
}

void Voronoi::geominit()
{
    double sn;
        update_k();


    freeinit(&efl, sizeof(Edge));
        update_k();

    nvertices = 0;
    nedges = 0;
    sn = (double)nsites + 4;
        update_k();

    sqrt_nsites = (int)sqrt(sn);
        update_k();

    deltay = ymax - ymin;
    deltax = xmax - xmin;
        update_k();

}

Halfedge*  Voronoi::HEcreate(Edge *e, int pm)
{
    Halfedge *answer;
        update_k();

    answer = (Halfedge *)getfree(&hfl);
        update_k();

    answer->ELedge = e;
        update_k();

    answer->ELpm = pm;
        update_k();

    answer->PQnext = (Halfedge *)NULL;
        update_k();

    answer->vertex = (Site *)NULL;
        update_k();

    answer->ELrefcnt = 0;
        update_k();

    return answer;
}

void  Voronoi::ELinsert(Halfedge *lb, Halfedge *newHe)
{
    newHe->ELleft = lb;
        update_k();

    newHe->ELright = lb->ELright;
        update_k();

    (lb->ELright)->ELleft = newHe;
        update_k();

    lb->ELright = newHe;
}

Halfedge *  Voronoi::ELgethash(int b)
{
    Halfedge *he;    update_k();


    if (b<0 || b >= ELhashsize)
        return((Halfedge *)NULL);
    he = ELhash[b];
    if (he == (Halfedge *)NULL || he->ELedge != (Edge *)-2)
        return (he);
    update_k();

    /* Hash table points to deleted half edge.  Patch as necessary. */
    ELhash[b] = (Halfedge *)NULL;
    if ((he->ELrefcnt -= 1) == 0)
        makefree((Freenode*)he, &hfl);
        update_k();

    return ((Halfedge *)NULL);
}

Halfedge *  Voronoi::ELleftbnd(VoronoiPoint *p)
{
    int i, bucket;
    Halfedge *he;
          update_k();

    bucket = (int)((p->x - xmin) / deltax * ELhashsize);
    update_k();

    if (bucket<0) bucket = 0;
    if (bucket >= ELhashsize) bucket = ELhashsize - 1;

    he = ELgethash(bucket);    update_k();

    if (he == (Halfedge *)NULL)
    {
        for (i = 1; 1; i += 1)
        {
            if ((he = ELgethash(bucket - i)) != (Halfedge *)NULL)
                break;
                update_k();

            if ((he = ELgethash(bucket + i)) != (Halfedge *)NULL)
                break;
                update_k();

        };
        totalsearch += i;
    };
    ntry += 1;
    if (he == ELleftend || (he != ELrightend && right_of(he, p)))
    {
        do
        {
            he = he->ELright;
                update_k();

        } while (he != ELrightend && right_of(he, p));
        he = he->ELleft;
    }
    else
        do
        {
            he = he->ELleft;
                update_k();

        } while (he != ELleftend && !right_of(he, p));

        if (bucket > 0 && bucket <ELhashsize - 1)
        {
            if (ELhash[bucket] != (Halfedge *)NULL)
            {
                ELhash[bucket]->ELrefcnt -= 1;
                    update_k();

            }
            ELhash[bucket] = he;
                update_k();

            ELhash[bucket]->ELrefcnt += 1;
                update_k();

        };
        return (he);
}
void Voronoi::ELdelete(Halfedge *he)
{
    (he->ELleft)->ELright = he->ELright;
    (he->ELright)->ELleft = he->ELleft;
        update_k();

    he->ELedge = (Edge *)-2;
        update_k();

    update_k();

}


Halfedge *  Voronoi::ELright(Halfedge *he)
{
        update_k();

    return (he->ELright);
}

Halfedge * Voronoi::ELleft(Halfedge *he)
{
        update_k();

    return (he->ELleft);
}


Site *  Voronoi::leftreg(Halfedge *he)
{
    if (he->ELedge == (Edge *)NULL)
        return(bottomsite);
        update_k();

    return(he->ELpm == 0 ? he->ELedge->Sites[0] : he->ELedge->Sites[1]);
}

Site * Voronoi::rightreg(Halfedge *he)
{
    if (he->ELedge == (Edge *)NULL)
        return(bottomsite);
        update_k();

    return(he->ELpm == 0 ? he->ELedge->Sites[1] : he->ELedge->Sites[0]);
}

Edge * Voronoi::bisect(Site *s1, Site *s2)
{
    double dx, dy, adx, ady;
    Edge *newedge;
        update_k();


    newedge = (Edge *)getfree(&efl);
        update_k();


    newedge->Sites[0] = s1; //store the sites that this edge is bisecting
    newedge->Sites[1] = s2;
        update_k();

    ref(s1);
        update_k();

    ref(s2);
        update_k();

    newedge->Vertices[0] = (Site *)NULL; //to begin with, there are no endpoints on the bisector - it goes to infinity
    newedge->Vertices[1] = (Site *)NULL;
         update_k();

    dx = s2->coord.x - s1->coord.x;         //get the difference in x dist between the sites
        update_k();

    dy = s2->coord.y - s1->coord.y;
            update_k();

    adx = dx>0 ? dx : -dx;                  //make sure that the difference in positive
        update_k();

    ady = dy>0 ? dy : -dy;
        update_k();

    newedge->c = (double)(s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5);//get the slope of the line
             update_k();

    if (adx>ady)
    {
        newedge->a = 1.0; newedge->b = dy / dx; newedge->c /= dx;//set formula of line, with x fixed to 1
    }
    else
    {
        newedge->b = 1.0; newedge->a = dx / dy; newedge->c /= dy;//set formula of line, with y fixed to 1
    };

    newedge->edgenbr = nedges;
        update_k();


    nedges += 1;
    return(newedge);
}

Site * Voronoi::intersect(Halfedge *el1, Halfedge *el2, VoronoiPoint *p)
{
    Edge *e1, *e2, *e;
        update_k();

    Halfedge *el;
        update_k();

    double d, xint, yint;
        update_k();

    int right_of_site;
        update_k();

    Site *v;

    e1 = el1->ELedge;
    e2 = el2->ELedge;
        update_k();

    if (e1 == (Edge*)NULL || e2 == (Edge*)NULL)
        return ((Site *)NULL);
        update_k();


    //if the two edges bisect the same parent, return null
    if (e1->Sites[1] == e2->Sites[1])
        return ((Site *)NULL);
        update_k();


    d = e1->a * e2->b - e1->b * e2->a;
    if (-1.0e-10<d && d<1.0e-10)
        return ((Site *)NULL);
        update_k();


    xint = (e1->c*e2->b - e2->c*e1->b) / d;
    yint = (e2->c*e1->a - e1->c*e2->a) / d;
        update_k();


    if ((e1->Sites[1]->coord.y < e2->Sites[1]->coord.y) ||
        (e1->Sites[1]->coord.y == e2->Sites[1]->coord.y &&
            e1->Sites[1]->coord.x < e2->Sites[1]->coord.x))
    {
        el = el1;
            update_k();

        e = e1;
            update_k();

    }
    else
    {
        el = el2;
            update_k();

        e = e2;
            update_k();

    };

    right_of_site = xint >= e->Sites[1]->coord.x;
    if ((right_of_site && el->ELpm == 0) || (!right_of_site && el->ELpm == 1))
        return ((Site *)NULL);
    update_k();

    //create a new site at the point of intersection - this is a new vector event waiting to happen
    v = (Site *)getfree(&sfl);
    v->refcnt = 0;
        update_k();

    v->coord.x = xint;
        update_k();

    v->coord.y = yint;
        update_k();

    return(v);
}
int Voronoi::right_of(Halfedge *el, VoronoiPoint *p)
{
    Edge *e;
    Site *topsite;
    int right_of_site, above, fast;
    double dxp, dyp, dxs, t1, t2, t3, yl;

    e = el->ELedge;
    topsite = e->Sites[1];
    right_of_site = p->x > topsite->coord.x;
    if (right_of_site && el->ELpm == 0) return(1);
    if (!right_of_site && el->ELpm == 1) return (0);

    if (e->a == 1.0)
    {
        dyp = p->y - topsite->coord.y;
        dxp = p->x - topsite->coord.x;
        fast = 0;
        if ((!right_of_site & (e->b<0.0)) | (right_of_site & (e->b >= 0.0)))
        {
            above = dyp >= e->b*dxp;
            fast = above;
        }
        else
        {
            above = p->x + p->y*e->b > e->c;
            if (e->b<0.0) above = !above;
            if (!above) fast = 1;
        };
        if (!fast)
        {
            dxs = topsite->coord.x - (e->Sites[0])->coord.x;
            above = e->b * (dxp*dxp - dyp*dyp) <
                dxs*dyp*(1.0 + 2.0*dxp / dxs + e->b*e->b);
            if (e->b<0.0) above = !above;
        };
    }
    else  /*e->b==1.0 */
    {
        yl = e->c - e->a*p->x;
        t1 = p->y - yl;
        t2 = p->x - topsite->coord.x;
        t3 = yl - topsite->coord.y;
        above = t1*t1 > t2*t2 + t3*t3;
    };
    return (el->ELpm == 0 ? above : !above);
}

void Voronoi::endpoint(Edge *e, int lr, Site * s)
{
    e->Vertices[lr] = s;
    ref(s);
    if (e->Vertices[1 - lr] == (Site *)NULL)
        return;
        update_k();


    clip_line(e);
            update_k();

    deref(e->Sites[0]);
    deref(e->Sites[1]);
        update_k();

    makefree((Freenode*)e, &efl);
        update_k();

}

double Voronoi::dist(Site *s, Site *t)
{
    double dx, dy;
        update_k();

    dx = s->coord.x - t->coord.x;
        update_k();

    dy = s->coord.y - t->coord.y;
        update_k();

    return (double)(sqrt(dx*dx + dy*dy));
}


void Voronoi::makevertex(Site *v)
{
    v->sitenbr = nvertices;
    nvertices += 1;
}


void Voronoi::deref(Site *v)
{
    v->refcnt -= 1;
        update_k();

    if (v->refcnt == 0)
        makefree((Freenode*)v, &sfl);
}

void  Voronoi::ref(Site *v)
{
    v->refcnt += 1;
        update_k();

}

void Voronoi::PQinsert(Halfedge *he, Site * v, double offset)
{
    Halfedge *last, *next;

    he->vertex = v;
        update_k();

    ref(v);
        update_k();

    he->ystar = (double)(v->coord.y + offset);
    last = &PQhash[PQbucket(he)];
        update_k();

    while ((next = last->PQnext) != (Halfedge *)NULL &&
        (he->ystar  > next->ystar ||
            (he->ystar == next->ystar && v->coord.x > next->vertex->coord.x)))
    {
        last = next;
    };
    he->PQnext = last->PQnext;
        update_k();

    last->PQnext = he;
        update_k();

    PQcount += 1;
}

void Voronoi::PQdelete(Halfedge *he)
{
    Halfedge *last;

    if (he->vertex != (Site *)NULL)
    {
        last = &PQhash[PQbucket(he)];
        while (last->PQnext != he)
            last = last->PQnext;

        last->PQnext = he->PQnext;
        PQcount -= 1;
        deref(he->vertex);
        he->vertex = (Site *)NULL;
    };
}

int  Voronoi::PQbucket(Halfedge *he)
{
    int bucket;

    bucket = (int)((he->ystar - ymin) / deltay * PQhashsize);
    if (bucket<0) bucket = 0;
    if (bucket >= PQhashsize) bucket = PQhashsize - 1;
    if (bucket < PQmin) PQmin = bucket;
    return(bucket);
}



int Voronoi::PQempty()
{
    return(PQcount == 0);
}


VoronoiPoint Voronoi::PQ_min()
{
    VoronoiPoint answer;

    while (PQhash[PQmin].PQnext == (Halfedge *)NULL) { PQmin += 1; };
    answer.x = PQhash[PQmin].PQnext->vertex->coord.x;
    answer.y = PQhash[PQmin].PQnext->ystar;
    return (answer);
}

Halfedge * Voronoi::PQextractmin()
{
    Halfedge *curr;

    curr = PQhash[PQmin].PQnext;
    PQhash[PQmin].PQnext = curr->PQnext;
    PQcount -= 1;
    return(curr);
}


bool Voronoi::PQinitialize()
{
    int i;

    PQcount = 0;
    PQmin = 0;
    PQhashsize = 4 * sqrt_nsites;
    PQhash = (Halfedge *)myalloc(PQhashsize * sizeof *PQhash);

    if (PQhash == 0)
        return false;

    for (i = 0; i<PQhashsize; i += 1) PQhash[i].PQnext = (Halfedge *)NULL;

    return true;
}


void Voronoi::freeinit(Freelist *fl, int size)
{
    fl->head = (Freenode *)NULL;
    fl->nodesize = size;
}

char *  Voronoi::getfree(Freelist *fl)
{
    int i;
    Freenode *t;

    if (fl->head == (Freenode *)NULL)
    {
        t = (Freenode *)myalloc(sqrt_nsites * fl->nodesize);

        if (t == 0)
            return 0;

        currentMemoryBlock->next = new FreeNodeArrayList;
        currentMemoryBlock = currentMemoryBlock->next;
        currentMemoryBlock->memory = t;
        currentMemoryBlock->next = 0;

        for (i = 0; i<sqrt_nsites; i += 1)
            makefree((Freenode *)((char *)t + i*fl->nodesize), fl);
    };
    t = fl->head;
    fl->head = (fl->head)->nextfree;
    return((char *)t);
}



void Voronoi::makefree(Freenode *curr, Freelist *fl)
{
    curr->nextfree = fl->head;
    fl->head = curr;
}

void  Voronoi::cleanup()
{
    if (sites != 0)
    {
        free(sites);
        sites = 0;
    }

    FreeNodeArrayList* current = 0, *prev = 0;

    current = prev = allMemoryList;

    while (current->next != 0)
    {
        prev = current;
        current = current->next;
        free(prev->memory);
        delete prev;
        prev = 0;
    }

    if (current != 0 && current->memory != 0)
    {
        free(current->memory);
        delete current;
    }

    allMemoryList = new FreeNodeArrayList;
    allMemoryList->next = 0;
    allMemoryList->memory = 0;
    currentMemoryBlock = allMemoryList;
}

void Voronoi::cleanupEdges()
{
    GraphEdge* geCurrent = 0, *gePrev = 0;
    geCurrent = gePrev = allEdges;

    while (geCurrent != 0 && geCurrent->next != 0)
    {
        gePrev = geCurrent;
        geCurrent = geCurrent->next;
        delete gePrev;
    }
    allEdges = 0;
}

void Voronoi::pushGraphEdge(double x1, double y1, double x2, double y2)
{
    GraphEdge* newEdge = new GraphEdge;
    newEdge->next = allEdges;
    allEdges = newEdge;
    newEdge->x1 = x1;
    newEdge->y1 = y1;
    newEdge->x2 = x2;
    newEdge->y2 = y2;
}


char * Voronoi::myalloc(unsigned n)
{
    char *t = 0;
    t = (char*)malloc(n);
    total_alloc += n;
    return(t);
}

void  Voronoi::line(double x1, double y1, double x2, double y2)
{
    pushGraphEdge(x1, y1, x2, y2);

}
void  Voronoi::clip_line(Edge *e)
{
    Site *s1, *s2;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0, temp = 0;;

    x1 = e->Sites[0]->coord.x;
    x2 = e->Sites[1]->coord.x;
    y1 = e->Sites[0]->coord.y;
    y2 = e->Sites[1]->coord.y;

    //if the distance between the two points this line was created from is less than 
    //the square root of 2, then ignore it
    if (sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites)
    {
        return;
    }
    pxmin = borderMinX;
    pxmax = borderMaxX;
    pymin = borderMinY;
    pymax = borderMaxY;

    if (e->a == 1.0 && e->b >= 0.0)
    {
        s1 = e->Vertices[1];
        s2 = e->Vertices[0];
    }
    else
    {
        s1 = e->Vertices[0];
        s2 = e->Vertices[1];
    };

    if (e->a == 1.0)
    {
        y1 = pymin;
        if (s1 != (Site *)NULL && s1->coord.y > pymin)
            y1 = s1->coord.y;
        if (y1>pymax)
            y1 = pymax;
        x1 = e->c - e->b * y1;
        y2 = pymax;
        if (s2 != (Site *)NULL && s2->coord.y < pymax)
            y2 = s2->coord.y;

        if (y2<pymin)
            y2 = pymin;
        x2 = (e->c) - (e->b) * y2;
        if (((x1> pxmax) & (x2>pxmax)) | ((x1<pxmin)&(x2<pxmin)))
            return;
        if (x1> pxmax)
        {
            x1 = pxmax; y1 = (e->c - x1) / e->b;
        };
        if (x1<pxmin)
        {
            x1 = pxmin; y1 = (e->c - x1) / e->b;
        };
        if (x2>pxmax)
        {
            x2 = pxmax; y2 = (e->c - x2) / e->b;
        };
        if (x2<pxmin)
        {
            x2 = pxmin; y2 = (e->c - x2) / e->b;
        };
    }
    else
    {
        x1 = pxmin;
        if (s1 != (Site *)NULL && s1->coord.x > pxmin)
            x1 = s1->coord.x;
        if (x1>pxmax)
            x1 = pxmax;
        y1 = e->c - e->a * x1;
        x2 = pxmax;
        if (s2 != (Site *)NULL && s2->coord.x < pxmax)
            x2 = s2->coord.x;
        if (x2<pxmin)
        {
            x2 = pxmin;
        }
        y2 = e->c - e->a * x2;
        if (((y1> pymax) & (y2>pymax)) | ((y1<pymin)&(y2<pymin)))
            return;
        if (y1> pymax)
        {
            y1 = pymax; x1 = (e->c - y1) / e->a;
        };
        if (y1<pymin)
        {
            y1 = pymin; x1 = (e->c - y1) / e->a;
        };
        if (y2>pymax)
        {
            y2 = pymax; x2 = (e->c - y2) / e->a;
        };
        if (y2<pymin)
        {
            y2 = pymin; x2 = (e->c - y2) / e->a;
        };
    };

    VEdge ee;
    ee.Left_Site = e->Sites[0]->coord;
    ee.Right_Site = e->Sites[1]->coord;
    ee.VertexA.x = x1;
    ee.VertexA.y = y1;
    ee.VertexB.x = x2;
    ee.VertexB.y = y2;

    total_edges.push_back(ee);
    line(x1, y1, x2, y2);
}

bool  Voronoi::voronoi(int triangulate)
{
    Site *newsite, *bot, *top, *temp, *p;
    Site *v;
    VoronoiPoint newintstar;
    int pm;
    Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
    Edge *e;

    PQinitialize();
    bottomsite = nextone();
    bool retval = ELinitialize();

    if (!retval)
        return false;

    newsite = nextone();
    while (1)
    {

        if (!PQempty())
            newintstar = PQ_min();
        if (newsite != (Site *)NULL && (PQempty() || newsite->coord.y < newintstar.y
            || (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x)))
        {
            lbnd = ELleftbnd(&(newsite->coord));
            rbnd = ELright(lbnd);
            bot = rightreg(lbnd);
            e = bisect(bot, newsite);
            bisector = HEcreate(e, 0);
            ELinsert(lbnd, bisector);

            if ((p = intersect(lbnd, bisector)) != (Site *)NULL)
            {
                PQdelete(lbnd);
                PQinsert(lbnd, p, dist(p, newsite));
            };
            lbnd = bisector;
            bisector = HEcreate(e, 1);
            ELinsert(lbnd, bisector);

            if ((p = intersect(bisector, rbnd)) != (Site *)NULL)
            {
                PQinsert(bisector, p, dist(p, newsite));
            };
            newsite = nextone();
        }
        else if (!PQempty())
        {
            lbnd = PQextractmin();
            llbnd = ELleft(lbnd);
            rbnd = ELright(lbnd);
            rrbnd = ELright(rbnd);
            bot = leftreg(lbnd);
            top = rightreg(rbnd);
            v = lbnd->vertex;
            makevertex(v);
            endpoint(lbnd->ELedge, lbnd->ELpm, v);
            endpoint(rbnd->ELedge, rbnd->ELpm, v);
            ELdelete(lbnd);
            PQdelete(rbnd);
            ELdelete(rbnd);
            pm = 0;
            if (bot->coord.y > top->coord.y)
            {
                temp = bot;
                bot = top;
                top = temp;
                pm = 1;
            }
            e = bisect(bot, top);
            bisector = HEcreate(e, pm);
            ELinsert(llbnd, bisector);
            endpoint(e, 1 - pm, v);
            deref(v);
            if ((p = intersect(llbnd, bisector)) != (Site *)NULL)
            {
                PQdelete(llbnd);
                PQinsert(llbnd, p, dist(p, bot));
            };

            //if right HE and the new bisector don't intersect, then reinsert it 
            if ((p = intersect(bisector, rrbnd)) != (Site *)NULL)
            {
                PQinsert(bisector, p, dist(p, bot));
            };
        }
        else break;
    };

    for (lbnd = ELright(ELleftend); lbnd != ELrightend; lbnd = ELright(lbnd))
    {
        e = lbnd->ELedge;

        clip_line(e);
    };

    cleanup();
    return true;
}

/* return a single in-storage site */
Site * Voronoi::nextone()
{
    Site *s;
    if (siteidx < nsites)
    {
        s = &sites[siteidx];
        siteidx += 1;
        return(s);
    }
    else
        return((Site *)NULL);
}

Voronoi* vdg;
vector<VoronoiPoint*> ver;
vector<VEdge> edges;

int main(int argc, char **argv)
{
    int offset = 200;
    int factor = 200;

    int iou = 0;
    for (vector<VoronoiPoint*>::iterator i = ver.begin(); i != ver.end(); i++)
        delete((*i));
    ver.clear();
    edges.clear();

    int n, seed;
    cout << "Enter the number of random sites: ";
    cin >> n ;
    cout << "Enter seed: ";
    cin >> seed;

    cout << "\nNumber of random sites: "<< n << endl;
    cout << "Seed: " << seed << endl;

    srand(seed);

    for(int i=0; i<n; i++)
    {
        double temp1 = (rand()%2000)/1000.0;
        double temp2 = (rand()%2000)/1000.0;

        ver.push_back(new VoronoiPoint(temp1, temp2));
    }
    /*
    ver.clear();
    ver.push_back(new VoronoiPoint(1.5, 1.5));
    ver.push_back(new VoronoiPoint(0.75, 0.99));
    */
    vdg = new Voronoi();
    double minY = -10;
    double maxY = 10;
    edges = vdg->ComputeVoronoiGraph(ver, minY, maxY);
    delete vdg;

    map <pair<double, double>, double> mapping;

    for(vector<VoronoiPoint*>::iterator i = ver.begin(); i!=ver.end(); i++)
    {
        mapping[ std::make_pair((*i)->x, (*i)->y) ] = 1E+20;
    }

    for(vector<VEdge>::iterator j = edges.begin(); j!=edges.end(); j++)
    {
        double distance_left = 100;
        double distance_right = 100;
        if(j->VertexA.x == j->VertexB.x)
        {
            distance_left = fabs(j->VertexA.x - j->Left_Site.x);
            distance_right = fabs(j->Right_Site.x - j->VertexA.x);
        }
        else
        {
            double slope = (j->VertexA.y - j->VertexB.y)/(j->VertexA.x - j->VertexB.x);
            double intercept = j->VertexA.y -(slope*j->VertexA.x);

            distance_left = fabs(j->Left_Site.y - slope*j->Left_Site.x - intercept)/sqrt((slope*slope) + 1);
            distance_right = fabs(j->Right_Site.y - slope*j->Right_Site.x - intercept)/sqrt((slope*slope) + 1);
        }
        mapping[std::make_pair(j->Left_Site.x, j->Left_Site.y)] = min(mapping[std::make_pair(j->Left_Site.x, j->Left_Site.y)], distance_left);
        mapping[std::make_pair(j->Right_Site.x, j->Right_Site.y)] = min(mapping[std::make_pair(j->Right_Site.x, j->Right_Site.y)], distance_right);
    }


    string svg_text_begin="<svg height=\"1000\" width=\"1000\" xmlns=\"http://www.w3.org/2000/svg\">\n<rect width=\"1000\" height=\"1000\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />\n";
    string point_string = "";
    for(vector<VoronoiPoint*>::iterator i = ver.begin(); i!=ver.end(); i++)
    {

        cout << (*i)->x << " " << (*i)->y << endl;

        stringstream x_string;
        stringstream y_string;
        stringstream r_string;
        x_string << factor*((*i)->x) + offset;
        y_string << factor*((*i)->y) + offset;
        r_string << (mapping[std::make_pair((*i)->x, (*i)->y)])*factor;
        point_string =  point_string + "<circle cx=\""+ x_string.str() +"\" cy=\""+ y_string.str() +"\" r=\"1\" stroke=\"black\" stroke-width=\"3\"/>\n";
        point_string =  point_string + "<circle cx=\""+ x_string.str() +"\" cy=\""+ y_string.str() +"\" r=\""+r_string.str()+"\" stroke=\"#0066ff\" stroke-width=\"0\" fill=\"#0066ff\" fill-opacity=\"0.4\"/>\n";
    }

    string edge_string = "";
    for(vector<VEdge>::iterator j = edges.begin(); j!=edges.end(); j++)
    {
        stringstream p1_x;
        stringstream p1_y;
        stringstream p2_x;
        stringstream p2_y;

        p1_x << j->VertexA.x*factor + offset;
        p1_y << j->VertexA.y*factor + offset;
        p2_x << j->VertexB.x*factor + offset;
        p2_y << j->VertexB.y*factor + offset;

        edge_string = edge_string + "<line x1=\""+p1_x.str()+"\" y1=\""+p1_y.str()+"\" x2=\""+p2_x.str()+"\" y2=\""+p2_y.str()+"\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
    }

    string svg_text_end = "</svg>";

    stringstream result;
    result << svg_text_begin << point_string << edge_string << svg_text_end << endl;

    cout << result.str();

    ofstream file;
    file.open("diagram.svg");
    file << result.str();
    file.close();

}