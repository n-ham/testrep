/*******************************************************************************

    Copyright (C) 2018 Nicholas C. Ham

    This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
    International License. See http://creativecommons.org/licenses/by-sa/4.0/

    For a discussion about this code and the latest version see
    https://gitlab.com/n-ham-paper-files/lattice-path-algorithms

*******************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <queue>
#include <vector>
#include <set>
#include <sstream>

#define PI 3.14159265

typedef std::pair<int, int> Point;

//addition operation for points
Point operator+(const Point &p1, const Point &p2)
{
    return Point(p1.first + p2.first, p1.second + p2.second);
}

//subtraction operation for points
Point operator-(const Point &p1, const Point &p2)
{
    return Point(p1.first - p2.first, p1.second - p2.second);
}

double distance(const Point &p1, const Point &p2)
{
    return sqrt(pow(p2.first-p1.first,2) + pow(p2.second-p1.second,2));
}

struct Step
{
    Point step_amount;
    int bend_amount;
    std::string bend_direction, colour;
};

//less than relation for steps
bool operator<(const Step& s1, const Step& s2)
{
    if(s1.step_amount < s2.step_amount)
        return 1;
    else if(s1.step_amount > s2.step_amount)
        return 0;

    if(s1.bend_amount < s2.bend_amount)
        return 1;
    else if(s1.bend_amount > s2.bend_amount)
        return 0;

    //must be equal
    return 0;
}

struct StepSet
{
    int no_steps;
    std::set<Step> steps;
};

/*
    input function for step sets

    input format:
        noSteps
        x-coord y-coord stepDirection bendAmount edgeColour

    sample input:
        3
        0 1 left 0 red
        1 0 left 0 blue
        2 2 left 30 green
*/
std::istream& operator>>(std::istream& is, StepSet& ss)
{
    Step inStep;

    is >> ss.no_steps;
    for(int s=0; s<ss.no_steps; s++)
    {
        is >> inStep.step_amount.first >> inStep.step_amount.second >> inStep.bend_direction >> inStep.bend_amount >> inStep.colour;
        ss.steps.insert(inStep);
    }

    return is;
}

//output function for step sets
std::ostream& operator<<(std::ostream& os, const StepSet& ss)
{
    os << "{";
    int i=0;
    for(auto step=ss.steps.begin(); step!=ss.steps.end(); step++, i++)
    {
        os << "(" << step->step_amount.first << "," << step->step_amount.second << ")";
        if(i+1 < ss.no_steps)
            os << ", ";
    }
    os << "}";

    return os;
}

struct Constraint
{
    std::string orderCondition;
    double slope, intercept;
};

struct ConstraintVector
{
    int noConstraints;
    std::vector<Constraint> constraints;
};

/*
    input function for constraint vectors

    input format:
        noConstraints
        constraintType slope intercept

    sample input:
        2
        <= 1 0
        > 0 0
*/
std::istream& operator>>(std::istream& is, ConstraintVector& cv)
{
    Constraint inConstraint;

    is >> cv.noConstraints;
    for(int s=0; s<cv.noConstraints; s++)
    {
        is >> inConstraint.orderCondition >> inConstraint.slope >> inConstraint.intercept;
        cv.constraints.push_back(inConstraint);
    }

    return is;
}

std::ostream& operator<<(std::ostream& os, const ConstraintVector& cv)
{
    os << "constraints:" << std::endl;
    int i=0;
    for(auto constraint=cv.constraints.begin(); constraint!=cv.constraints.end(); constraint++, i++)
        os << "y " << constraint->orderCondition << " " << constraint->slope << "x + " << constraint->intercept << std::endl;

    return os;
}

typedef std::set<Point> EndPointSet;

struct Edge
{
    Point source, target;
    int bendAmount;
    std::string bendDirection, colour;

    Edge(const Point& Source, const Point& Target, int BendAmount, const std::string& BendDirection, const std::string& Colour)
    {
        source = Source;
        target = Target;
        bendAmount = BendAmount;
        bendDirection = BendDirection;
        colour = Colour;
    }
};

//less than relation on edges
bool operator<(const Edge& e1, const Edge& e2)
{
    if(e1.source < e2.source)
        return 1;
    else if(e1.source > e2.source)
        return 0;

    if(e1.target < e2.target)
        return 1;
    else if(e1.target > e2.target)
        return 0;

    //must be equal
    return 0;
}

struct Graph
{
    EndPointSet vertices;
    std::map<Point, long long int> labels;
    std::set<Edge> edges;
};

//output function for graphs
//NOTE: can't input const Graph& and use graph.labels[*vertex]
std::ostream& operator<<(std::ostream& os, Graph& graph)
{
    os << "    \\begin{tikzpicture}" << std::endl;
	os << "        \\tikzstyle{vertex}=[circle, draw=black, fill=white, inner sep = 0.06cm]" << std::endl;
	os << "        \\tikzstyle{distvertex}=[circle, draw=black, fill=lightgray, inner sep = 0.06cm]" << std::endl;

    for(auto vertex=graph.vertices.begin(); vertex != graph.vertices.end(); vertex++)
    {
        os << "        \\node";
        if(vertex->first == 0 && vertex->second == 0)
            os << "[distvertex] ";
        else
            os << "[vertex] ";
        os << "(" << vertex->first << "-" << vertex->second << ") at (" << vertex->first << ", " << vertex->second << ") {$" << graph.labels[*vertex] << "$};" << std::endl;

    }

    for(auto edge=graph.edges.begin(); edge != graph.edges.end(); edge++)
    {
        /*
            if you get 'Dimension too large' error(s), try not outputting
            edges between vertices with large labels, especially vertices
            close together.
        */
        /*if(graph.labels[edge->source] > 99 && graph.labels[edge->target] > 99 && distance(edge->source, edge->target) < 2.0)
            continue;*/
        /*if((graph.labels[edge->source] > 999 || graph.labels[edge->target] > 999) && distance(edge->source, edge->target) < 2.0)
            continue;*/
        os << "        \\draw [->-=0.55, ultra thick, " << edge->colour << "] (" << edge->source.first << "-" << edge->source.second << ") to [bend " << edge->bendDirection << "=" << edge->bendAmount << "] (" << edge->target.first << "-" << edge->target.second << ");" << std::endl;
    }

	os << "    \\end{tikzpicture}" << std::endl;

    return os;
}


/*
    algorithm 1

    inputs: 1. finite step set X \subseteq \mathbb{Z}^2_{\times};
            2. finite vector of constraints C;
            3. q \in \mathbb{N}; and
            4. bool value indicating whether to save TikZ output to output.txt.

    output: 1. the end point sets \mathscr{A}_X^C(i) for i = 0, 1, ..., q; and
            2. the graphs \Lambda_X^C(i) for i = 0, 1, ..., q.

    NOTE: output is put in to the last two input variables endPointSets and graphs
*/
void algorithm1(const StepSet& stepSet, const ConstraintVector& constraintVector, int q, const bool &saveTikZ, std::vector<EndPointSet> &endPointSets, std::vector<Graph> &graphs)
{
    endPointSets = std::vector<EndPointSet>(q+1, EndPointSet());
    graphs = std::vector<Graph>(q+1, Graph());

    //handles i=0 case
    endPointSets[0].insert(Point(0, 0));
    graphs[0].vertices.insert(Point(0, 0));
    graphs[0].labels[Point(0, 0)] = 0;

    //i = 1, ..., q cases
    for(int i=1; i<=q; i++)
    {
        graphs[i].vertices = graphs[i-1].vertices;
        graphs[i].edges = graphs[i-1].edges;

        for(auto point=endPointSets[i-1].begin(); point != endPointSets[i-1].end(); point++)
        {
            for(auto step=stepSet.steps.begin(); step != stepSet.steps.end(); step++)
            {
                bool valid = 1;
                Point nPoint = *point + step->step_amount;

                //this is a good place to hard-code constraints on x, eg. x >= 0 or x <= 10
                /*if(nPoint.first < 0) //x >= 0
                    valid = 0;*/

                for(auto constraint=constraintVector.constraints.begin(); constraint != constraintVector.constraints.end(); constraint++)
                {
                    if(constraint->orderCondition == "<=" && nPoint.second > constraint->slope*nPoint.first + constraint->intercept)
                        valid = 0;
                    else if(constraint->orderCondition == "<" && nPoint.second >= constraint->slope*nPoint.first + constraint->intercept)
                        valid = 0;
                    else if(constraint->orderCondition == ">=" && nPoint.second < constraint->slope*nPoint.first + constraint->intercept)
                        valid = 0;
                    else if(constraint->orderCondition == ">" && nPoint.second <= constraint->slope*nPoint.first + constraint->intercept)
                        valid = 0;
                }

                if(valid)
                {
                    endPointSets[i].insert(*point + step->step_amount);
                    graphs[i].vertices.insert(*point + step->step_amount);
                    graphs[i].edges.insert(Edge(*point, *point + step->step_amount, step->bend_amount, step->bend_direction, step->colour));
                }
            }
        }

        for(auto point=graphs[i].vertices.begin(); point != graphs[i].vertices.end(); point++)
        {
            if(endPointSets[i].count(*point))
                graphs[i].labels[*point] = i;
            else
                graphs[i].labels[*point] = graphs[i-1].labels[*point];
        }
    }

    if(saveTikZ)
    {
        std::ofstream ofs("./output.txt");
        for(int i=0; i<=q; i++)
            ofs << graphs[i] << std::endl << "    \\vspace{7mm}" << std::endl << std::endl;
        ofs.close();
    }
}


/*
    algorithm 2

    inputs: 1. finite step set X \subseteq \mathbb{Z}^2_{\times}; and
            2. bool value indicating whether to save TikZ output to output.txt.

    output: 1. whether X satisfies FPP or IPP; and
            2. vector/Point u such that the line through O and perpendicular
               to u witnesses the LC, with u pointing towards the half-plane
               containing x.

    NOTE: output is put in to the last two input variables property and u
*/
void algorithm2(const StepSet& stepSet, const bool &saveTikZ, std::string &property, Point &u)
{
    std::priority_queue<std::pair<double, Point>, std::vector<std::pair<double, Point> >, std::greater<std::pair<double, Point> > > prique;
    std::vector<std::pair<double, Point> > orderedSteps;
    std::set<double> angles;

    property = "";

    //calculates the angle \alpha_i from the positive x-axis to OA_i
    for(auto step=stepSet.steps.begin(); step != stepSet.steps.end(); step++)
    {
        if(step->step_amount.first > 0 && step->step_amount.second == 0) //positive x-axis
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is on pos x-axis with angle 0" << std::endl;
            prique.push(std::pair<double, Point>(0.0, step->step_amount));
            angles.insert(0.0);
        }
        else if(step->step_amount.first > 0 && step->step_amount.second > 0) //first quadrant
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is in first quadrant with angle " << atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI << std::endl;
            prique.push(std::pair<double, Point>(atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI, step->step_amount));
            angles.insert(atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI);
        }
        else if(step->step_amount.first == 0 && step->step_amount.second > 0) //positive y-axis
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is on pos y-axis with angle 90" << std::endl;
            prique.push(std::pair<double, Point>(90.0, step->step_amount));
            angles.insert(90.0);
        }
        else if(step->step_amount.first < 0 && step->step_amount.second > 0) //second quadrant
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is in second quadrant with angle " << atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI << std::endl;
            prique.push(std::pair<double, Point>(90.0 + atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI, step->step_amount));
            angles.insert(90.0 + atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI);
        }
        else if(step->step_amount.first < 0 && step->step_amount.second == 0) //negative x-axis
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is on neg x-axis with angle 180" << std::endl;
            prique.push(std::pair<double, Point>(180.0, step->step_amount));
            angles.insert(180.0);
        }
        else if(step->step_amount.first < 0 && step->step_amount.second < 0) //third quadrant
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is in third quadrant with angle " << atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI << std::endl;
            prique.push(std::pair<double, Point>(180.0 + atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI, step->step_amount));
            angles.insert(180.0 + atan((double)step->step_amount.second/(double)step->step_amount.first)*180/PI);
        }
        else if(step->step_amount.first == 0 && step->step_amount.second < 0) //negative y-axis
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is on neg y-axis with angle 270" << std::endl;
            prique.push(std::pair<double, Point>(270.0, step->step_amount));
            angles.insert(270.0);
        }
        else if(step->step_amount.first > 0 && step->step_amount.second < 0) //fourth quadrant
        {
            //std::cout << "step (" << step->step_amount.first << ", " << step->step_amount.second << ") is in fourth quadrant with angle " << atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI << std::endl;
            prique.push(std::pair<double, Point>(270.0 + atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI, step->step_amount));
            angles.insert(270.0 + atan((double)-step->step_amount.first/(double)step->step_amount.second)*180/PI);
        }
    }

    //std::cout << std::endl;
    while(!prique.empty())
    {
        //std::cout << "(" << prique.top().second.first << ", " << prique.top().second.second << ") " << prique.top().first << std::endl;
        orderedSteps.push_back(prique.top());
        prique.pop();
    }

    if(angles.size() == 1) //\alpha_{j_1} = ... = \alpha_{j_n}
    {
        property = "FPP";
        u = Point(orderedSteps[0].second.first, orderedSteps[0].second.second);
    }
    else
    {
        double beta_sum = 0;
        for(size_t s=0; s<orderedSteps.size(); s++)
        {
            double beta;
            if(s+1 < orderedSteps.size())
                beta = orderedSteps[s+1].first - orderedSteps[s].first;
            else
                beta = 360.0 + orderedSteps[0].first - orderedSteps[s].first;
            //std::cout << "beta_" << s << " = " << beta << std::endl;

            if(beta > 180.0)
            {
                property = "FPP";

                if(s+1 < orderedSteps.size())
                    u = Point(orderedSteps[s].second.second - orderedSteps[s+1].second.second,
                              orderedSteps[s+1].second.first - orderedSteps[s].second.first);
                else
                    u = Point(orderedSteps[s].second.second - orderedSteps[0].second.second,
                              orderedSteps[0].second.first - orderedSteps[s].second.first);

                break;
            }

            beta_sum += beta;
            //std::cout << "sum beta_i = " << beta_sum << std::endl;

            if(beta_sum >= 180.0)
            {
                property = "IPP";
                break;
            }
        }

        if(property == "")
            property = "IPP";
    }

    if(saveTikZ)
    {
        std::ofstream ofs("output.txt");

        ofs << "\\begin{tikzpicture}" << std::endl;

        ofs << "    \\fill[blue!20] (-4,4)--(-4,-4)--(4,-4)--(4,4);" << std::endl;

        int i=1;
        for(auto step=orderedSteps.begin(); step != orderedSteps.end(); step++, i++)
        {
            ofs << "    \\fill (" <<  step->second.first << "," << step->second.second << ") circle( 0.15);" << std::endl;
            ofs << "    \\draw[ultra thick] (0,0)--(" << step->second.first  << "," << step->second.second << ");" << std::endl;
            if(step->second.first >= 0 && step->second.second >= 0)
                ofs << "    \\node () at (" <<  step->second.first + 0.3 << "," <<  step->second.second + 0.3 << ") {$A_" <<  i << "$};" << std::endl;
            else if(step->second.first < 0 && step->second.second >= 0)
                ofs << "    \\node () at (" <<  step->second.first - 0.3 << "," <<  step->second.second + 0.3 << ") {$A_" <<  i << "$};" << std::endl;
            else if(step->second.first < 0 && step->second.second < 0)
                ofs << "    \\node () at (" <<  step->second.first - 0.3 << "," <<  step->second.second - 0.3 << ") {$A_" <<  i << "$};" << std::endl;
            else if(step->second.first >= 0 && step->second.second < 0)
                ofs << "    \\node () at (" <<  step->second.first + 0.3 << "," <<  step->second.second - 0.3 << ") {$A_" <<  i << "$};" << std::endl;
        }
        ofs << std::endl;

        ofs << "    \\draw[ultra thick, red, -{latex}] (0,0)--(" << u.first  << "," << u.second << ");" << std::endl;
        if(u.first >= 0 && u.second >= 0)
            ofs << "    \\node[red] () at (" <<  u.first + 0.2 << "," <<  u.second + 0.2 << ") {$\\textbf{u}$};" << std::endl;
        else if(u.first < 0 && u.second >= 0)
            ofs << "    \\node[red] () at (" <<  u.first - 0.2 << "," <<  u.second + 0.2 << ") {$\\textbf{u}$};" << std::endl;
        else if(u.first < 0 && u.second < 0)
            ofs << "    \\node[red] () at (" <<  u.first - 0.2 << "," <<  u.second - 0.2 << ") {$\\textbf{u}$};" << std::endl;
        else if(u.first >= 0 && u.second < 0)
            ofs << "    \\node[red] () at (" <<  u.first + 0.2 << "," <<  u.second - 0.2 << ") {$\\textbf{u}$};" << std::endl;
        ofs << std::endl;

        if(u.second == 0)
        {
            ofs << "    \\draw[ultra thick, blue, <->] (0,4)--(0,-4);" << std::endl;
            ofs << "    \\node[blue] () at (0.4,-3.6) {$\\mathscr{L}$};" << std::endl;
        }
        else
        {
            double slope = -(double)u.first/(double)u.second;
            if((slope >= 0 && slope <= 1) || (slope < 0 && slope >= -1))
            {
                ofs << "    \\draw[ultra thick, blue, <->] (-4," << -4.0*slope << ")--(4," << 4.0*slope << ");" << std::endl;
                ofs << "    \\node[blue] () at (3.7," << 3.7*slope + 0.4 << ") {$\\mathscr{L}$};" << std::endl;
            }
            else
            {
                ofs << "    \\draw[ultra thick, blue, <->] (" << -4.0/slope << ",-4)--(" << 4.0/slope << ",4);" << std::endl;
                ofs << "    \\node[blue] () at (" << 3.7/slope + 0.4 << ",3.7) {$\\mathscr{L}$};" << std::endl;
            }
        }
        ofs << std::endl;

        ofs << "    \\draw[ultra thick, fill=white] (0,0) circle (0.15);" << std::endl;
        ofs << "    \\node (0-0) at (-.3,-.3) {$O$};" << std::endl;

        ofs << "\\end{tikzpicture}" << std::endl;

        ofs.close();
    }
}

/*
    algorithm 3

    inputs: 1. finite step set X \subseteq \mathbb{Z}^2_{\times};
            2. finite vector of constraints;
            3. q \in \mathbb{N}; and
            4. bool value indicating whether to save TikZ output to output.txt.

    output: graph with
            1. vertex set \mathscr{A}_X[<=q];
            2. vertices labelled with numbers \pi_X(A) for each A \in \mathscr{A}_X[<=q]; and
            3. edges A->A+B for all A \in \mathscr{A}_X[<q] and B \in X.

    NOTE: output is put in to the last input variable graph
*/
void algorithm3(const StepSet& stepSet, const ConstraintVector& constraintVector, int q, const bool &saveTikZ, Graph &graph)
{
    double mu1 = 500.0, mu2 = -500.0;
    int Q;

    std::string property;
    Point u;

    algorithm2(stepSet, 0, property, u);

    for(auto step=stepSet.steps.begin(); step != stepSet.steps.end(); step++)
    {
        mu1 = std::min(mu1, (double)u.first*step->step_amount.first + u.second*step->step_amount.second);
        mu2 = std::max(mu2, (double)u.first*step->step_amount.first + u.second*step->step_amount.second);
    }

    Q = floor(mu2/mu1*q);

    std::vector<EndPointSet> endPointSets;
    std::vector<Graph> graphs;
    algorithm1(stepSet, constraintVector, Q, 0, endPointSets, graphs);

    std::vector<EndPointSet> endPointSets2(q+1, EndPointSet());

    for(auto point=graphs[q].vertices.begin(); point != graphs[q].vertices.end(); point++)
    {
        if(graphs[Q].labels[*point] <= q)
        {
            endPointSets2[graphs[Q].labels[*point] ].insert(*point);
            graph.vertices.insert(*point);
        }
    }

    graph.labels[Point(0, 0)] = 1;
    for(int i=1; i<=q; i++)
    {
        for(auto point=endPointSets2[i].begin(); point != endPointSets2[i].end(); point++)
        {
            graph.labels[*point] = 0;

            for(auto step=stepSet.steps.begin(); step != stepSet.steps.end(); step++)
            {
                if(graph.vertices.count(*point - step->step_amount))
                {
                    graph.labels[*point] += graph.labels[*point - step->step_amount];
                    graph.edges.insert(Edge(*point - step->step_amount, *point, step->bend_amount, step->bend_direction, step->colour));
                }
            }
        }
    }

    if(saveTikZ)
    {
        std::ofstream ofs("./output.txt");
        ofs << graph << std::endl;
        ofs.close();
    }
}
