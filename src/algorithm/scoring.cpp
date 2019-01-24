#include <sstream>
#include <vector>
#include <fstream>
#include <set>
#include <nlohmann/json.hpp>
#include <map>

#include "algorithm/electron_optimisation.hpp"
#include "utils/options.hpp"

template <class T, class ICT = std::vector<T>> struct CartesianProduct {
  using type = T;
  using innerType = ICT;
  
  using innerIter = typename innerType::const_iterator;
  
  struct innerIters {
    innerIter first, second, third;
    innerIters(innerIter a, innerIter b, innerIter c) : first(a), second(b), third(c) { }
  };
  
  using innerItersC = std::vector<innerIters>;
  
  innerItersC iters;
  bool finished;
  
  CartesianProduct() = delete;
  
  template <class outerIter>
  CartesianProduct(outerIter begin, outerIter end) : finished(false) {
    for (; begin != end; ++begin)
      iters.emplace_back(begin->begin(), begin->end(), begin->begin());
    if (iters.empty())
      finished = true;
  }
  
  CartesianProduct(innerType &a, innerType &b) : finished(false) {
    iters.emplace_back(a.begin(), a.end(), a.begin());
    iters.emplace_back(b.begin(), b.end(), b.begin());
    if (a.begin() == a.end() || b.begin() == b.end())
      finished = true;
  }
  
  // Returns true if c contains a product, false otherwise
  bool operator()(innerType &c) {
    c.clear();
    if (finished)
      return false;
    c.reserve(iters.size());
    for (auto it : iters)
      c.push_back(*(it.third));
    
    for (auto it = iters.begin();;) {
      ++(it->third);
      if (it->third == it->second) {
        if (it + 1 == iters.end())
          finished = true;
        else {
          it->third = it->first;
          ++it;
        }
      } else
        break;
    }
    return true;
  }
};

namespace indigox {
  
  Scoring ElectronOpt::scores_ = Scoring();
  String ElectronOpt::loaded_file = "";
  
  using json = nlohmann::json;
  using VString = std::vector<String>;
  // Checks if all the required keys are present, and if any unsupported keys are present
  String KeyChecker(VString& required, VString& optional, json& data) {
    std::set<String> req(required.begin(), required.end()), opt(optional.begin(), optional.end());
    std::stringstream error_stream;
    std::unordered_map<String, json> key_value_pairs;
    data.get_to(key_value_pairs);
    
    bool added_additional = false;
    
    for (auto& key_value : key_value_pairs) {
      auto req_pos = req.find(key_value.first);
      auto opt_pos = opt.find(key_value.first);
      if (req_pos != req.end()) req.erase(req_pos);
      else if (opt_pos == opt.end() && !key_value.second.is_null()) {
        if (!added_additional) {
          error_stream << "Additional keys=";
          added_additional = true;
        }
        error_stream << key_value.first << ",";
      }
    }
    if (!req.empty()) {
      error_stream << "Missing required keys=";
      for (String key : req) error_stream << key << ",";
    }
    return error_stream.str();
  }
  
  bool value_test(json& data, String type) {
    if (type == "str") return data.is_string();
    else if (type == "uint") return data.is_number_unsigned();
    else if (type == "int") return data.is_number_integer();
    else if (type == "obj") return data.is_object();
    else if (type == "array") return data.is_array();
    else throw std::runtime_error("Invalid value test type");
  }
  
  // Checks if all values have correct type
  String ValueChecker(VString& req_keys, VString& req_type, VString& opt_keys, VString& opt_type, json& data) {
    std::stringstream error_stream;
    bool added_header = false;
    
    for (size_t i = 0; i < req_keys.size(); ++i) {
      if (!value_test(data[req_keys[i]], req_type[i])) {
        if (!added_header) {
          error_stream << "Invalid types: ";
          added_header = true;
        }
        error_stream << "(key=" << req_keys[i] << ",expected=" << req_type[i] << "),";
      }
    }
    
    for (size_t i = 0; i < opt_keys.size(); ++i) {
      if (!data[opt_keys[i]].is_null() && !value_test(data[opt_keys[i]], opt_type[i])) {
        if (!added_header) {
          error_stream << "Invalid types: ";
          added_header = true;
        }
        error_stream << "(key=" << req_keys[i] << ",expected=" << req_type[i] << "),";
      }
    }
    return error_stream.str();
  }
  
  void AtomMatchingType::SanityCheck(json &data) {
    VString required_keys;
    VString required_type;
    VString optional_keys = {"element", "degree", "fc"};
    VString optional_type = {"str", "uint", "int"};
    
    // Check is object
    if (!data.is_object()) {
      throw std::runtime_error("Require JSON object type");
    }
    
    // Check keys and value types
    String errors = KeyChecker(required_keys, optional_keys, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    errors = ValueChecker(required_keys, required_type, optional_keys, optional_type, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    
    // Check sensibility of provided keys
    int count = 0;
    for (String& key : optional_keys) {
      if (!data[key].is_null()) ++count;
    }
    if (count == 0) {
      std::stringstream error;
      error << "No data provided: " << data;
      throw std::runtime_error(error.str());
    }
  }
  
  void AtomScoring::SanityCheck(json& data) {
    VString required_keys = {"atom"};
    VString required_type = {"obj"};
    VString optional_keys = {"neighbours", "scores", "score"};
    VString optional_type = {"array", "array", "uint"};
    
    // Check is object
    if (!data.is_object()) {
      throw std::runtime_error("Require JSON object type");
    }
    
    // Check keys and value types
    String errors = KeyChecker(required_keys, optional_keys, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    errors = ValueChecker(required_keys, required_type, optional_keys, optional_type, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    
    // Check sensibility of provided values, scores need valence and score values
    if (!data["scores"].is_null() && !data["score"].is_null()) {
      throw std::runtime_error("Only one scoring type supported");
    }
    if (data["scores"].is_null() && data["score"].is_null()) {
      throw std::runtime_error("No scores provided");
    }
    if (!data["scores"].is_null()) {
      VString req_score_key = {"valence", "score"};
      VString req_score_typ = {"uint", "uint"};
      VString empty;
      int score_count = 0;
      for (json& score : data["scores"]) {
        score_count++;
        errors = KeyChecker(req_score_key, empty, score);
        if (!errors.empty()) {
          throw std::runtime_error(errors);
        }
        errors = ValueChecker(req_score_key, req_score_typ, empty, empty, score);
        if (!errors.empty()) {
          throw std::runtime_error(errors);
        }
      }
      
      if (!score_count) {
        throw std::runtime_error("No scores provided");
      }
    }
  }
  
  void BondScoring::SanityCheck(json& data) {
    VString required_keys = {"atomA", "atomB"};
    VString required_type = {"obj", "obj"};
    VString optional_keys = {"neighboursA", "neighboursB", "scores", "score"};
    VString optional_type = {"array", "array", "array", "uint"};
    
    // Check is object
    if (!data.is_object()) {
      throw std::runtime_error("Require JSON object type");
    }
    
    // Check keys and value types
    String errors = KeyChecker(required_keys, optional_keys, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    errors = ValueChecker(required_keys, required_type, optional_keys, optional_type, data);
    if (!errors.empty()) {
      throw std::runtime_error(errors);
    }
    
    // Check sensibility of provided values, scores need valence and score values
    if (!data["scores"].is_null() && !data["score"].is_null()) {
      throw std::runtime_error("Only one scoring type supported");
    }
    if (data["scores"].is_null() && data["score"].is_null()) {
      throw std::runtime_error("No scores provided");
    }
    if (!data["scores"].is_null()) {
      VString req_score_key = {"order", "score"};
      VString req_score_typ = {"uint", "uint"};
      VString empty;
      int score_count = 0;
      for (json& score : data["scores"]) {
        score_count++;
        errors = KeyChecker(req_score_key, empty, score);
        if (!errors.empty()) {
          throw std::runtime_error(errors);
        }
        errors = ValueChecker(req_score_key, req_score_typ, empty, empty, score);
        if (!errors.empty()) {
          throw std::runtime_error(errors);
        }
      }
      
      if (!score_count) {
        throw std::runtime_error("No scores provided");
      }
    }
  }
  
#define ATOMFC_POSITIVE 15
#define ATOMFC_NEGATIVE 14
  AtomMatchingType::AtomMatchingType(json& data) : element(0), degree(0), formal_charge(0){
    SanityCheck(data);
    json& test_data = data["element"];
    if (!test_data.is_null()) {
      PeriodicTable_p pt = PeriodicTable::GetInstance();
      String symb = test_data.get<String>();
      if (symb == "X") {
        element.set();
      } else {
        uint8_t Z = pt->GetElement(symb)->GetAtomicNumber();
        element.set(Z);
      }
    } else {
      element.set();
    }
    
    test_data = data["degree"];
    if (!test_data.is_null()) {
      degree.set(test_data.get<size_t>());
    } else {
      degree.set();
    }
    
    test_data = data["fc"];
    if (!test_data.is_null()) {
      int8_t fc = test_data.get<int8_t>();
      formal_charge.set(abs(fc));
      if (fc < 0) formal_charge.set(ATOMFC_NEGATIVE);
      if (fc > 0) formal_charge.set(ATOMFC_POSITIVE);
    } else {
      formal_charge.set();
    }
  }
  
  bool AtomMatchingType::Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex &vertex) {
    ElnVertProp* p = graph->GetProperties(vertex);
    if (!element.test(p->atomic_number)) return false;
    if (!degree.test(graph->Degree(vertex))) return false;
    if (!formal_charge.test(abs(p->formal_charge))) return false;
    if (p->formal_charge < 0 && !formal_charge[ATOMFC_NEGATIVE]) return false;
    if (p->formal_charge > 0 && !formal_charge[ATOMFC_POSITIVE]) return false;
    return true;
  }
  
#define MAXIMUM_VALENCE 16
  AtomScoring::AtomScoring(json& data) {
    SanityCheck(data);
    
    // required infos
    atom = AtomMatchingType(data["atom"]);
    
    // optional infos
    json& test_data = data["scores"];
    if (!test_data.is_null()) {
      for (json& score : test_data) {
        scoring[score["valence"].get<uint8_t>()] = score["score"].get<Score>();
      }
    }
    
    test_data = data["score"];
    if (!test_data.is_null()) {
      for (uint8_t valence = 0; valence < MAXIMUM_VALENCE; ++valence) {
        scoring[valence] = test_data.get<Score>();
      }
    }
    
    test_data = data["neighbours"];
    if (!test_data.is_null()) {
      for (json& atom : test_data) neighbours.emplace_back(atom);
    }
  }
  
  bool AtomScoring::Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex &vertex) {
    if (!atom.Matches(graph, vertex)) return false;
    if (neighbours.empty()) return true;
    
    ElnVertProp* prop = graph->GetProperties(vertex);
    auto nbrs = graph->GetNeighbours(vertex);
    std::vector<std::vector<size_t>> matching_positions;
    std::vector<ElnVertex> neighbour_vertices;
    
    for (; nbrs.first != nbrs.second; ++nbrs.first) {
      ElnVertProp* edge_prop = graph->GetProperties(*nbrs.first);
      if (edge_prop->id.first == prop->id.first) {
        neighbour_vertices.emplace_back(graph->GetVertex(std::make_pair(edge_prop->id.second, edge_prop->id.second)));
      } else {
        neighbour_vertices.emplace_back(graph->GetVertex(std::make_pair(edge_prop->id.first, edge_prop->id.first)));
      }
    }
    
    for (size_t test_nbr = 0; test_nbr < neighbours.size(); ++test_nbr) {
      matching_positions.emplace_back(std::vector<size_t>());
      for (size_t vert_nbr = 0; vert_nbr < neighbour_vertices.size(); ++vert_nbr) {
        if (neighbours[test_nbr].Matches(graph, neighbour_vertices[vert_nbr])) {
          matching_positions[test_nbr].push_back(vert_nbr);
        }
      }
      if (matching_positions[test_nbr].empty()) return false;
    }
    
    CartesianProduct<size_t> product(matching_positions.begin(), matching_positions.end());
    std::vector<size_t> test_product;
    while (product(test_product)) {
      std::set<size_t> size_test(test_product.begin(), test_product.end());
      if (size_test.size() == test_product.size()) return true;
    }
    
    return false;
  }
 
#define MAXIMUM_BOND_ORDER 16
  BondScoring::BondScoring(json& data) {
    SanityCheck(data);
    
    // required infos
    atom_a = AtomMatchingType(data["atomA"]);
    atom_b = AtomMatchingType(data["atomB"]);
    
    // optional infos
    json& test_data = data["scores"];
    if (!test_data.is_null()) {
      for (json& score : test_data) {
        scoring[score["order"].get<uint8_t>()] = score["score"].get<Score>();
      }
    }
    
    test_data = data["score"];
    if (!test_data.is_null()) {
      for (uint8_t order = 0; order < MAXIMUM_BOND_ORDER; ++order) {
        scoring[order] = test_data.get<Score>();
      }
    }
    
    test_data = data["neighboursA"];
    if (!test_data.is_null()) {
      for (json& atom : test_data) neighbours_a.emplace_back(atom);
    }
    test_data = data["neighboursB"];
    if (!test_data.is_null()) {
      for (json& atom : test_data) neighbours_b.emplace_back(atom);
    }
  }
  
  bool BondScoring::Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex &vertex) {
    /// \todo Make matching test the neighbours as well
    ElnVertProp* prop = graph->GetProperties(vertex);
    ElnVertex source = graph->GetVertex(std::make_pair(prop->id.first, prop->id.first));
    ElnVertex target = graph->GetVertex(std::make_pair(prop->id.second, prop->id.second));
    if (atom_a.Matches(graph, source) && atom_b.Matches(graph, target)) return true;
    std::swap(source, target);
    return (atom_a.Matches(graph, source) && atom_b.Matches(graph, target));
  }
  
  void Scoring::LoadScoreFile(String path) {
    atom_scores.clear();
    bond_scores.clear();
    potential_atom_matches.clear();
    potential_bond_matches.clear();
    json j;
    std::ifstream input(path);
    input >> j;
    
    if (!j.is_array()) throw std::runtime_error("Input JSON file must contain an array of entries");
    
    for (auto& entry : j) {
      if (!entry.is_object()) throw std::runtime_error("Each entry must be a JSON object.");
      
      if (!entry["atom"].is_null()) atom_scores.emplace_back(entry);
      else if (!entry["atomA"].is_null() && !entry["atomB"].is_null()) bond_scores.emplace_back(entry);
      else {
        std::stringstream error;
        error << "Unable to determine entry type: " << entry;
        throw std::runtime_error(error.str());
      }
    }
    
    std::vector<uint8_t> supported_elements = {1,6,7,8,9,15,16,17,35};
    for (uint8_t element : supported_elements) {
      for (size_t i = 0; i < atom_scores.size(); ++i) {
        if (atom_scores[i].atom.element.test(element)) {
          potential_atom_matches[element].push_back(i);
        }
      }
      for (size_t i = 0; i < bond_scores.size(); ++i) {
        if (bond_scores[i].atom_a.element.test(element) || bond_scores[i].atom_b.element.test(element)) {
          potential_bond_matches[element].push_back(i);
        }
      }
    }
  }
  
  Score Scoring::GetScore(std::shared_ptr<ElectronGraph> graph, ElnVertex &vertex) {
    ElnVertProp* prop = graph->GetProperties(vertex);
    if (prop->id.first == prop->id.second) {
      AtomScoring* matcher = nullptr;
      if (potential_atom_matches[prop->atomic_number].empty()) {
        for (AtomScoring& test : atom_scores) {
          if (test.Matches(graph, vertex)) {
            matcher = &test;
            break;
          }
        }
      } else {
        for (size_t i : potential_atom_matches[prop->atomic_number]) {
          if (atom_scores[i].Matches(graph, vertex)) {
            matcher = &atom_scores[i];
            break;
          }
        }
      }
      
      if (matcher == nullptr) return Options::AssignElectrons::INF;
      uint8_t valence = 0;
      for (auto nbrs = graph->GetNeighbours(vertex); nbrs.first != nbrs.second; ++nbrs.first) {
        ElnVertProp* p = graph->GetProperties(*nbrs.first);
        valence += p->electron_count + p->pre_placed;
      }
      valence /= 2;
      
      auto pos = matcher->scoring.find(valence);
      return (pos == matcher->scoring.end()) ? Options::AssignElectrons::INF : pos->second;
    } else {
//      return 0;
      BondScoring* matcher = nullptr;
      ElnVertex a = graph->GetVertex(std::make_pair(prop->id.first, prop->id.first));
      ElnVertex b = graph->GetVertex(std::make_pair(prop->id.second, prop->id.second));
      
      uint8_t Za = graph->GetProperties(a)->atomic_number;
      uint8_t Zb = graph->GetProperties(b)->atomic_number;

      uint8_t working;
      if (potential_bond_matches[Za].size() <= potential_bond_matches[Zb].size() && !potential_bond_matches[Za].empty()) {
        working = Za;
      } else if (potential_bond_matches[Zb].size() < potential_bond_matches[Za].size() && !potential_bond_matches[Zb].empty()) {
        working = Zb;
      } else {
        working = 0;
      }
      
      
      if (working == 0) {
        for (BondScoring& test : bond_scores) {
          if (test.Matches(graph, vertex)) {
            matcher = &test;
            break;
          }
        }
      } else {
        for (size_t i : potential_bond_matches[working]) {
          if (bond_scores[i].Matches(graph, vertex)) {
            matcher = &bond_scores[i];
            break;
          }
        }
      }
      
      if (matcher == nullptr) return Options::AssignElectrons::INF;
      
      auto pos = matcher->scoring.find((prop->electron_count + prop->pre_placed)/2);
      
      return (pos == matcher->scoring.end()) ? Options::AssignElectrons::INF : pos->second;
    }
  }
  
}
