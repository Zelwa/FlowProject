/*

LOISEL Lucas - 21711075
GUERET Morgan - 21700094

Soit notre graphe G = (V,E)

1 - Le nombre de supporteurs empruntant le bus reliant la ville «i» avec la ville «j» ne doit
pas excéder le nombre maximum de places disponibles pour les supporteurs.


Pour tout arc (i, j) de E : 0<= X(i, j) <= CAP(i, j)
CAP(i,j) : représentant le nombre maximum de places dans le bus en provenance de la ville i vers la ville j.


-----------

2 - Le nombre des supporteurs arrivant à une ville intermédiaire est égale au nombre de
supporteurs repartant à partir de cette ville.


D'après la conservation des flots :
Pour la ville 3 on a : X(2, 3) = X(3, 4) + X(3, D)
c'est-à-dire pour la ville 3 : sum(X(i, 3)) = sum(X(3, j))

Donc pour chaque ville v on a : sum(X(i, v)) = sum(X(v, j)) 

-----------

3 - L’expression linéaire de la fonction objectif à maximiser.


max : Z = sum(X(i, D))
D : la ville de destination


*/

#include "ortools/linear_solver/linear_solver.h"
#include "ortools/graph/max_flow.h"

#include <iostream>
#include <fstream>
#include <string>

typedef long int FlowQuantity;
typedef int NodeIndex;
typedef std::vector<std::pair<std::pair<NodeIndex, NodeIndex>, FlowQuantity>> vectorARC;

/**
* Meta Programme permettant de délimiter une string 
* à l'aide d'un délimitateur (utilisé dans un parser)
*
* @parm[in] str - ligne traité  
* @parm[in] cont - conteneur des caractères après delimitation
* @parm[in] delim - caractère de délimitation
**/
template <class Container>
void split2(const std::string& str, Container& cont, char delim = ' ') {
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    } // while
} // split2

// Bibliothèque operations_research contenant le solver
namespace operations_research {
  
  /**
  * Modélisation d'un problème linéaire de maximisation de flot
  * tiré d'ORTOOLS et des connaissances du cours.
  * 
  * @parm[in] nb_node - le nombre de noeuds de notre flot
  * @parm[in] arcs - le vecteur des arcs/capacités de notre flot
  **/
  void app_flow(int nb_node, vectorARC arcs) {

    const int num_nodes = nb_node;
    StarGraph graph(num_nodes, arcs.size());
    MaxFlow typed_flow(&graph, 0, num_nodes - 1);
    for (const auto& it : arcs) {
      ArcIndex arc = graph.AddArc(it.first.first, it.first.second);
      typed_flow.SetArcCapacity(arc, it.second);
    } // for

    LOG(INFO) << "Solving max flow with: " << graph.num_nodes() << " nodes, and "
              << graph.num_arcs() << " arcs.";

    typed_flow.Solve();
    if (MaxFlow::OPTIMAL != typed_flow.status()) {
      LOG(FATAL) << "Solving the max flow is not optimal!";
    } // if
    FlowQuantity total_flow = typed_flow.GetOptimalFlow();
    LOG(INFO) << "Maximum flow: " << total_flow;
    LOG(INFO) << "";
    LOG(INFO) << " Arc  : Flow / Capacity";
    for (int i = 0; i < arcs.size(); ++i) {
      LOG(INFO) << graph.Tail(i) << " -> " << graph.Head(i) << ": "
                << typed_flow.Flow(i) << " / " << typed_flow.Capacity(i);
    } // for
  } // app_flow

  /**
  * Modélisation d'un problème linéaire de maximisation de flot
  * 
  **/
  void app(int nb_arcs, 
            int source, int destination, 
            std::vector<int> begin_node, 
            std::vector<int> end_node, 
            std::vector<int> cap) {
    MPSolver solver("simple_lp_program", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    double infini = MPSolver::infinity();

    // matrice des arcs
    int arcs[nb_arcs][2] = { 0 };

    // on numérote les arcs en partant de 0
    if (source > 0){
      for(int i = 0; i < nb_arcs; i++){
          arcs[i][0] = begin_node[i] - source;
          arcs[i][1] = end_node[i] - source;
      }
      destination = destination - source;
      source = source - source;
    }
    else{
      for(int i = 0; i < nb_arcs; i++){
          arcs[i][0] = begin_node[i];
          arcs[i][1] = end_node[i];
      }
    }

    // Variables : matrice des flots
    MPVariable *flow[nb_arcs] = { 0 };
    for (int i = 0; i < nb_arcs; i++) {
      std::stringstream nom;
      nom << "x_" << (i);
      flow[i] = solver.MakeIntVar(0.0, infini, nom.str());
    }

    // Contraintes de capacité
    for (int i = 0; i < nb_arcs; i++) {
      std::stringstream nom;
      nom << "C_transport-" << (i + 1);

      LinearExpr tmp = flow[i];
      solver.MakeRowConstraint(tmp <= cap[i]);
    }


    // Contraintes inflows == outflows
    for(int i = 0; i < nb_arcs; i++) {

      LinearExpr inflows;
      for (int k = 0; k < nb_arcs; k++) {
        if (arcs[k][1] == i and arcs[k][1] != destination) {
           inflows += flow[k];
        }
      }

      LinearExpr outflows;
      for (int k = 0; k < nb_arcs; k++) {
        if (arcs[k][0] == i and arcs[k][0] != source) {
           outflows += flow[k];
        }
      }

      solver.MakeRowConstraint(inflows == outflows);
    }

    LOG(INFO) << "Number of variables = " << solver.NumVariables();
    LOG(INFO) << "Number of constraints = " << solver.NumConstraints() << std::endl;
    LOG(INFO) << "";

    // Fonction Objectif
    MPObjective *const objective = solver.MutableObjective();
    for (int k = 0; k < nb_arcs; k++) {
      if (arcs[k][0] == source) {
        objective->SetCoefficient(flow[k], 1);
      }
    } 
    objective->SetMaximization();

    solver.Solve();

    LOG(INFO) << "***** Solution *****" << std::endl;
    LOG(INFO) << "Valeur optimale de F = " << objective->Value() << std::endl;

    for (int i = 0; i < nb_arcs; ++i) {
      LOG(INFO) << "Arc : " << arcs[i][0] << " -> " << arcs[i][1]
                << " / flux = " << flow[i]->solution_value() 
                << " / " << cap[i];
    }
  } // app

} // operations_research

int main (int argc, char* argv[]) {
  using namespace std;
  cout << endl;

  string line;
  if (argc > 1){
    LOG(INFO) << "FILENAME : " << argv[1] << endl;
    ifstream file ("graphs/" + string(argv[1]));

    string lp_type;
    int nb_node, nb_arcs, source, destination;
    
    vectorARC arcs;
    std::vector<int> begin_node;
    std::vector<int> end_node;
    std::vector<int> cap;

    if (file.is_open())
    {
      while ( getline (file,line) )
      {
        std::vector<std::string> elems;
        const char delim = ' ';

        switch(line[0]){
          case 'c':
            break;
          case 'p':
            split2(line, elems, delim);
            lp_type = elems[1];
            nb_node = stoi(elems[2]);
            nb_arcs = stoi(elems[3]);
            break;
          case 'n':
            split2(line, elems, delim);
            if (elems[2] == "s"){
              source = stoi(elems[1]);
            }
            if (elems[2] == "t"){
              destination = stoi(elems[1]);
            }
            break;
          case 'a':
            split2(line, elems, delim);
            std::pair<NodeIndex, NodeIndex> p1;
            if (source > 0) {
              p1.first = stoi(elems[1]) - source;
              p1.second = stoi(elems[2]) - source;
            }
            else{
              p1.first = stoi(elems[1]);
              p1.second = stoi(elems[2]);
            }
            
            std::pair<std::pair<NodeIndex, NodeIndex>, FlowQuantity> contain;
            contain.first = p1;
            contain.second = stoi(elems[3]);
            arcs.push_back(contain);

            begin_node.push_back(stoi(elems[1]));
            end_node.push_back(stoi(elems[2]));
            cap.push_back(stoi(elems[3]));
            break;
        } // switch 
      } // while
      file.close();
    } // if
    else { cout << "Unable to open file" << endl; return EXIT_FAILURE;} // else
    LOG(INFO) << "";
    LOG(INFO) << "Fonction utilisant l'api : ortools/graph/max_flow.h" << std::endl;
    LOG(INFO) << "";
    operations_research::app_flow(nb_node, arcs);
    LOG(INFO) << "";
    LOG(INFO) << "#########################";
    LOG(INFO) << "";
    LOG(INFO) << "Fonction utilisant l'api : ortools/linear_solver/linear_solver.h" << std::endl;
    LOG(INFO) << "";
    operations_research::app(nb_arcs, source, destination, begin_node, end_node, cap);
    LOG(INFO) << "";
    LOG(INFO) << "END" << std::endl;
  }
  else{
    cout << "Filename not found, missing 1 argument" << endl; return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
} // main
