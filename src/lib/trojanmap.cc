#include "trojanmap.h"

//-----------------------------------------------------
// TODO: Students should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id. If id does not exist, return
 * -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(const std::string &id) { 
  return 0;
}

/**
 * GetLon: Get the longitude of a Node given its id. If id does not exist,
 * return -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(const std::string &id) {
  return 0;
}

/**
 * GetName: Get the name of a Node given its id. If id does not exist, return
 * "NULL".
 *
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(const std::string &id) {
  return "";
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node. If id does not exist, return
 * an empty vector.
 *
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(const std::string &id) {
  return {};
}

/**
 * GetID: Given a location name, return the id.
 * If the node does not exist, return an empty string.
 * The location name must be unique, which means there is only one node with the name.
 *
 * @param  {std::string} name          : location name
 * @return {std::string}               : id
 */
std::string TrojanMap::GetID(const std::string &name) {
  std::string res = "";
  for (const auto &pair : data) {
    if (pair.second.name == name) {
      return pair.first;  // Found matching name, return the node's ID
    }
  }
  return res; 
}
  

/**
 * GetPosition: Given a location name, return the position. If id does not
 * exist, return (-1, -1).
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);
  if (name.empty()){
    return results;
  }
  for (const auto& pair : data){
    if (pair.second.name == name){
      return {pair.second.lat, pair.second.lon};
    }
  }
  
  return results;
}

/**
 * CalculateEditDistance: Calculate edit distance between two location names
 * @param  {std::string} a          : first string
 * @param  {std::string} b          : second string
 * @return {int}                    : edit distance between two strings
 */
int TrojanMap::CalculateEditDistance(std::string a, std::string b) {
  int a_len = a.size();
  int b_len = b.size();
  
  std::vector<std::vector<int>> dynamic_programming(a_len+1, std::vector<int>(b_len+1, 0));

  std::transform(a.begin(), a.end(), a.begin(), ::tolower);
  std::transform(b.begin(), b.end(), b.begin(), ::tolower);

  for (int i=0; i<=a_len; i++) dynamic_programming[i][0] = i;
  for (int j=0; j<=b_len; j++) dynamic_programming[0][j] = j;

  for (int i=1; i<=a_len; i++){
    for (int j=1; j<=b_len; j++){
      if( a[i-1] == b[j-1]){
        dynamic_programming[i][j] = dynamic_programming[i-1][j-1];
      }
      else{
        dynamic_programming[i][j] = std::min({
          dynamic_programming[i-1][j-1],
          dynamic_programming[i][j-1],
          dynamic_programming[i-1][j]
        })+1;
      }
    }
  }
  return dynamic_programming[a_len][b_len];
}

/**
 * FindClosestName: Given a location name, return the name with the smallest edit
 * distance.
 *
 * @param  {std::string} name          : location name
 * @return {std::string} tmp           : the closest name
 */
std::string TrojanMap::FindClosestName(std::string name) {
  std::string tmp = ""; // Start with a dummy word
  int min_distance = INT_MAX;

  for (const auto &pair : data) {
    std::string pair_name = pair.second.name;
    
    if (pair_name.empty()) continue;

    int distance = CalculateEditDistance(name, pair_name);

    if (distance < min_distance) {
      min_distance = distance;
      tmp = pair_name;
    }
  }
  return tmp;
}

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix. The function should be case-insensitive.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
  std::vector<std::string> results;
  if (name.empty()){
    return results;
  }
  // convert input to lowercase
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);

  for (const auto& pair : data){
    std::string location_name = pair.second.name;
    std::string location_name_lower = location_name;
    std::transform(location_name_lower.begin(), location_name_lower.end(), location_name_lower.begin(), ::tolower);
  
    if(location_name_lower.find(name)==0){
      results.push_back(location_name); // return original name
    }
  }
  return results;
}

/**
 * GetAllCategories: Return all the possible unique location categories, i.e.
 * there should be no duplicates in the output.
 *
 * @return {std::vector<std::string>}  : all unique location categories
 */
std::vector<std::string> TrojanMap::GetAllCategories() {
  std::unordered_set<std::string> loc_categories;
  
  for (const auto& pair : data) {
    for (const auto& attribute : pair.second.attributes) {
      loc_categories.insert(attribute);
    }
  }
  std::vector<std::string> result(loc_categories.begin(), loc_categories.end());
  return result;
}

/**
 * GetAllLocationsFromCategory: Return all the locations of the input category (i.e.
 * 'attributes' in data.csv). If there is no location of that category, return
 * (-1, -1). The function should be case-insensitive.
 *
 * @param  {std::string} category         : category name (attribute)
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetAllLocationsFromCategory(
    std::string category) {
  std::vector<std::string> res;

  std::transform(category.begin(), category.end(), category.begin(), ::tolower);
  for (const auto& pair : data) {
    for (const auto& attribute : pair.second.attributes) {
      std::string attribute_lower = attribute;
      std::transform(attribute_lower.begin(), attribute_lower.end(), attribute_lower.begin(), ::tolower);

      if (attribute_lower == category) {
        res.push_back(pair.first);
        break;
      }
    }
  }
  return res;
}

/**
 * GetLocationRegex: Given the regular expression of a location's name, your
 * program should first check whether the regular expression is valid, and if so
 * it returns all locations that match that regular expression.
 *
 * @param  {std::regex} location name      : the regular expression of location
 * names
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetLocationRegex(std::regex location) {
  std::vector<std::string> results;

  for (const auto& pair : data) {
    const std::string& loc_name = pair.second.name;
    if (!loc_name.empty() && std::regex_match(loc_name, location)) {
      results.push_back(pair.first);
    }
  }
  return results;
}

/**
 * CalculateDistance: Get the distance between 2 nodes.
 * We have provided the code for you. Please do not need to change this function.
 * You can use this function to calculate the distance between 2 nodes.
 * The distance is in mile.
 * The distance is calculated using the Haversine formula.
 * https://en.wikipedia.org/wiki/Haversine_formula
 * 
 * @param  {std::string} a  : a_id
 * @param  {std::string} b  : b_id
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const std::string &a_id,
                                    const std::string &b_id) {
  // Do not change this function
  Node a = data[a_id];
  Node b = data[b_id];
  double dlon = (b.lon - a.lon) * M_PI / 180.0;
  double dlat = (b.lat - a.lat) * M_PI / 180.0;
  double p = pow(sin(dlat / 2), 2.0) + cos(a.lat * M_PI / 180.0) *
                                           cos(b.lat * M_PI / 180.0) *
                                           pow(sin(dlon / 2), 2.0);
  double c = 2 * asin(std::min(1.0, sqrt(p)));
  return c * 3961;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations
 * inside the vector.
 * We have provided the code for you. Please do not need to change this function.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  // Do not change this function
  double sum = 0;
  for (int i = 0; i < int(path.size()) - 1; i++) {
    sum += CalculateDistance(path[i], path[i + 1]);
  }
  return sum;
}

/**
 * CalculateShortestPath_Dijkstra: Given 2 locations, return the shortest path
 * which is a list of id. Hint: Use priority queue.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(
    std::string location1_name, std::string location2_name) {
  std::vector<std::string> path;

  std::string start_id = GetID(location1_name);
  std::string end_id = GetID(location2_name);
  if (start_id.empty() || end_id.empty()) return path;
  std::unordered_map<std::string, double> distance;
  std::unordered_map<std::string, std::string> previous;

  auto cmp = [&](const std::string &a, const std::string &b) {
    return distance[a] > distance[b];  // Using comparison function to find
  };                                   // smallest distance
  std::priority_queue<std::string, std::vector<std::string>, decltype(cmp)> pq(cmp);
  //Using priority queue to pick smallest distance

  for (auto &pair : data) {
    distance[pair.first] = DBL_MAX; //set all distance to infinity
  }

  distance[start_id] = 0; 
  pq.push(start_id);

  while (!pq.empty()) {
    std::string current = pq.top();
    pq.pop();

    if (current == end_id) break;

    for (const auto &neighbor_id : data[current].neighbors) {
      double d = distance[current] + CalculateDistance(current, neighbor_id);
      if (d < distance[neighbor_id]) {
        distance[neighbor_id] = d; //update shortest path
        previous[neighbor_id] = current;
        pq.push(neighbor_id);
      }
    }
  }

  std::string temp = end_id;
  while (temp != start_id && previous.count(temp)) {
    path.push_back(temp); //collect path from end -> start
    temp = previous[temp];
  }
  if (temp == start_id) {
    path.push_back(start_id);
    std::reverse(path.begin(), path.end()); //reverse to start -> end
  }
  else {
    path.clear();  // no path found
  }

  return path;
}

/**
 * CalculateShortestPath_Bellman_Ford: Given 2 locations, return the shortest
 * path which is a list of id. Hint: Do the early termination when there is no
 * change on distance.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(
    std::string location1_name, std::string location2_name) {
  std::vector<std::string> path;
  std::string start_id = GetID(location1_name);
  std::string end_id = GetID(location2_name);

  if (start_id.empty() || end_id.empty()) return path;

  std::unordered_map<std::string, double> distance;
  std::unordered_map<std::string, std::string> previous;

  for (auto &node : data) {
    distance[node.first] = DBL_MAX; // infinity
  }
  distance[start_id] = 0;

  int n = data.size(); 

  for (int i = 0; i < n - 1; i++) { 
    bool updated = false;

    for (auto &node : data) {
      std::string u = node.first; // Current node ID

      for (const std::string &v : node.second.neighbors) {
        double w = CalculateDistance(u, v); // Get weight of edge u -> v

        // Relax the edge only if the current known path is shorter
        if (distance[u] != DBL_MAX && distance[u] + w < distance[v]) {
          distance[v] = distance[u] + w;  // Update shortest distance to v
          previous[v] = u;
          updated = true;
        }
      }
    }
    if (!updated) break;
  }

  std::string temp = end_id;
  while (temp != start_id && previous.count(temp)) {
    path.push_back(temp);
    temp = previous[temp];
  }
  if (temp == start_id) {
    path.push_back(start_id);
    std::reverse(path.begin(), path.end());
  } else {
    path.clear();  // no path found
  }

  return path;
}

/**
 * Traveling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path, 
 *                                                                      for example: {10.3, {{0, 1, 2, 3, 4, 0}, {0, 1, 2, 3, 4, 0}, {0, 4, 3, 2, 1, 0}}},
 *                                                                      where 10.3 is the total distance, 
 *                                                                      and the first vector is the path from 0 and travse all the nodes and back to 0,
 *                                                                      and the second vector is the path shorter than the first one,
 *                                                                      and the last vector is the shortest path.
 */
// Please use brute force to implement this function, ie. find all the permutations.
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Brute_force(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

// Please use backtracking to implement this function
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Backtracking(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

// Hint: https://en.wikipedia.org/wiki/2-opt
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_2opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

// This is optional
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_3opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

/**
 * Given CSV filename, it read and parse locations data from CSV file,
 * and return locations vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_locations.csv"
 *   File content:
 *    Name
 *    Ralphs
 *    KFC
 *    Chick-fil-A
 *   Output: ['Ralphs', 'KFC', 'Chick-fil-A']
 * @param  {std::string} locations_filename     : locations_filename
 * @return {std::vector<std::string>}           : locations
 */
std::vector<std::string> TrojanMap::ReadLocationsFromCSVFile(
    std::string locations_filename) {
  std::vector<std::string> location_names_from_csv;
  std::fstream fin;
  fin.open(locations_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, word)) {
    location_names_from_csv.push_back(word);
  }
  fin.close();
  return location_names_from_csv;
}

/**
 * Given CSV filenames, it read and parse dependencise data from CSV file,
 * and return dependencies vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_dependencies.csv"
 *   File content:
 *     Source,Destination
 *     Ralphs,Chick-fil-A
 *     Ralphs,KFC
 *     Chick-fil-A,KFC
 *   Output: [['Ralphs', 'Chick-fil-A'], ['Ralphs', 'KFC'], ['Chick-fil-A', 'KFC']]
 * @param  {std::string} dependencies_filename     : dependencies_filename
 * @return {std::vector<std::vector<std::string>>} : dependencies
 */
std::vector<std::vector<std::string>> TrojanMap::ReadDependenciesFromCSVFile(
    std::string dependencies_filename) {
  std::vector<std::vector<std::string>> dependencies_from_csv;
  std::fstream fin;
  fin.open(dependencies_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);
    std::vector<std::string> dependency;
    while (getline(s, word, ',')) {
      dependency.push_back(word);
    }
    dependencies_from_csv.push_back(dependency);
  }
  fin.close();
  return dependencies_from_csv;
}

/**
 * DeliveringTrojan: Given a vector of location names, it should return a
 * sorting of nodes that satisfies the given dependencies. If there is no way to
 * do it, return a empty vector.
 *
 * @param  {std::vector<std::string>} locations                     : locations
 * @param  {std::vector<std::vector<std::string>>} dependencies     : prerequisites
 * @return {std::vector<std::string>} results                       : results
 */
std::vector<std::string> TrojanMap::DeliveringTrojan(
    std::vector<std::string> &locations,
    std::vector<std::vector<std::string>> &dependencies) {
  std::vector<std::string> result;
  std::unordered_map<std::string, std::vector<std::string>> list;
  std::unordered_map<std::string, int> incoming;

  for (const auto &loc : locations) {
    incoming[loc] = 0; //Initialize count
  }

  for (const auto &dep : dependencies) {
    std::string u = dep[0];  // comes first
    std::string v = dep[1];  // comes after
    list[u].push_back(v);
    incoming[v]++; // v have one incoming
  }

  std::queue<std::string> q;
  for (const auto &pair : incoming) {
    if (pair.second == 0) {
      q.push(pair.first); // starting node
    }
  }

  while (!q.empty()) {
    std::string current = q.front();
    q.pop();
    result.push_back(current);

    for (const auto &neighbor : list[current]) {
      incoming[neighbor]--;
      if (incoming[neighbor] == 0) {
        q.push(neighbor); // add to the queue until no prerequsite
      }
    }
  }

  if (result.size() != locations.size()) {
    result.clear();  // check if there is a cycle or invalid input
  }

  return result;     
}

/**
 * inSquare: Give a id retunr whether it is in square or not.
 *
 * @param  {std::string} id            : location id
 * @param  {std::vector<double>} square: four vertexes of the square area
 * @return {bool}                      : in square or not
 */
bool TrojanMap::inSquare(std::string id, std::vector<double> &square) {
  return true;
}


/**
 * GetSubgraph: Give four vertexes of the square area, return a list of location
 * ids in the squares
 *
 * @param  {std::vector<double>} square         : four vertexes of the square
 * area
 * @return {std::vector<std::string>} subgraph  : list of location ids in the
 * square
 */
std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
  // include all the nodes in subgraph
  std::vector<std::string> subgraph;
  for (const auto &pair : data) {
    const auto &node = pair.second;

    // Check if the node inside the square [0]:left [1]:right [2]:upper [3]:lower
    if (node.lon >= square[0] && node.lon <= square[1] &&
        node.lat >= square[3] && node.lat <= square[2]) {
      subgraph.push_back(node.id);
    }
  }
  return subgraph;
}

/**
 * Cycle Detection: Given four points of the square-shape subgraph, return true
 * if there is a cycle path inside the square, false otherwise.
 *
 * @param {std::vector<std::string>} subgraph: list of location ids in the
 * square
 * @param {std::vector<double>} square: four vertexes of the square area
 * @return {bool}: whether there is a cycle or not
 */
bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square) {
  std::unordered_set<std::string> node_set(subgraph.begin(), subgraph.end());
  std::unordered_map<std::string, bool> visited;
  std::function<bool(std::string, std::string)> dfs = [&](std::string current, std::string parent) {
    visited[current] = true;
    for (const auto &neighbor : data[current].neighbors) {
      if (node_set.count(neighbor)) { // Only coount neighbors in subgraph
        if (!visited[neighbor]) {
          if (dfs(neighbor, current)) return true;
        } else if (neighbor != parent) {
          // already visited but not parent
          return true;
        }
      }
    }
    return false;
  };

  // Check all unvisited nodes
  for (const auto &node_id : subgraph) {
    if (!visited[node_id]) {
      if (dfs(node_id, "")) return true;
    }
  }

  return false; // No cycles found
}

/**
 * FindNearby: Given a class name C, a location name L and a number r,
 * find all locations in class C on the map near L with the range of r and
 * return a vector of string ids
 *
 * @param {std::string} className: the name of the class
 * @param {std::string} locationName: the name of the location
 * @param {double} r: search radius
 * @param {int} k: search numbers
 * @return {std::vector<std::string>}: location name that meets the requirements
 */
std::vector<std::string> TrojanMap::FindNearby(std::string attributesName, std::string name, double r, int k) {
  std::vector<std::string> res;
  return res;
}

/**
 * Shortest Path to Visit All Nodes: Given a list of locations, return the shortest
 * path which visit all the places and no need to go back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::vector<std::string> }      : the shortest path
 */
std::vector<std::string> TrojanMap::TrojanPath(
      std::vector<std::string> &location_names) {
    std::vector<std::string> res;
    return res;
}

/**
 * Given a vector of queries, find whether there is a path between the two locations with the constraint of the gas tank.
 *
 * @param  {std::vector<std::pair<double, std::vector<std::string>>>} Q : a list of queries
 * @return {std::vector<bool> }      : existence of the path
 */
std::vector<bool> TrojanMap::Queries(const std::vector<std::pair<double, std::vector<std::string>>>& q) {
    std::vector<bool> ans(q.size());
    return ans;
}

/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * We have provided the code for you. Please do not need to change this function.
 */
void TrojanMap::CreateGraphFromCSVFile() {
  // Do not change this function
  std::fstream fin;
  fin.open("src/lib/data.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '{'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '}'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        if (isalpha(word[0])) n.attributes.insert(word);
        if (isdigit(word[0])) n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}
