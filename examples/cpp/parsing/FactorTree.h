// Copyright (c) 2012 Andre Martins
// All Rights Reserved.
//
// This file is part of AD3 2.1.
//
// AD3 2.1 is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// AD3 2.1 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with AD3 2.1.  If not, see <http://www.gnu.org/licenses/>.

#ifndef FACTOR_TREE_H_
#define FACTOR_TREE_H_

#include "ad3/GenericFactor.h"

namespace AD3 {

class Arc {
 public:
  Arc(int h, int m) : h_(h), m_(m) {}
  ~Arc() {}

  int head() { return h_; }
  int modifier() { return m_; }

 private:
  int h_;
  int m_;
};

class FactorTree : public GenericFactor {
 public:
  FactorTree() {}
  virtual ~FactorTree() { ClearActiveSet(); }

  int RunCLE(const std::vector<double>& scores,
             std::vector<int> *heads,
             double *value);

  // Compute the score of a given assignment.
  // Note: additional_log_potentials is empty and is ignored.
  void Maximize(const std::vector<double> &variable_log_potentials,
                const std::vector<double> &additional_log_potentials,
                Configuration &configuration,
                double *value) {
    std::vector<int>* heads = static_cast<std::vector<int>*>(configuration);
    RunCLE(variable_log_potentials, heads, value);
  }

  // Compute the score of a given assignment.
  // Note: additional_log_potentials is empty and is ignored.
  void Evaluate(const std::vector<double> &variable_log_potentials,
                const std::vector<double> &additional_log_potentials,
                const Configuration configuration,
                double *value) {
    const std::vector<int> *heads = static_cast<const std::vector<int>*>(configuration);
    // Heads belong to {0,1,2,...}
    *value = 0.0;
    for (int m = 1; m < heads->size(); ++m) {
      int h = (*heads)[m];
      int index = index_arcs_[h][m];
      *value += variable_log_potentials[index];
    }
  }

  // Given a configuration with a probability (weight),
  // increment the vectors of variable and additional posteriors.
  // Note: additional_log_potentials is empty and is ignored.
  void UpdateMarginalsFromConfiguration(
    const Configuration &configuration,
    double weight,
    std::vector<double> *variable_posteriors,
    std::vector<double> *additional_posteriors) {
    const std::vector<int> *heads = static_cast<const std::vector<int>*>(configuration);
    for (int m = 1; m < heads->size(); ++m) {
      int h = (*heads)[m];
      int index = index_arcs_[h][m];
      (*variable_posteriors)[index] += weight;
    }
  }

  // Count how many common values two configurations have.
  int CountCommonValues(const Configuration &configuration1,
                        const Configuration &configuration2) {
    const std::vector<int> *heads1 = static_cast<const std::vector<int>*>(configuration1);
    const std::vector<int> *heads2 = static_cast<const std::vector<int>*>(configuration2);
    int count = 0;
    for (int i = 1; i < heads1->size(); ++i) {
      if ((*heads1)[i] == (*heads2)[i]) {
        ++count;
      }
    }
    return count;
  }

  // Check if two configurations are the same.
  bool SameConfiguration(
    const Configuration &configuration1,
    const Configuration &configuration2) {
    const std::vector<int> *heads1 = static_cast<const std::vector<int>*>(configuration1);
    const std::vector<int> *heads2 = static_cast<const std::vector<int>*>(configuration2);
    for (int i = 1; i < heads1->size(); ++i) {
      if ((*heads1)[i] != (*heads2)[i]) return false;
    }
    return true;
  }

  // Delete configuration.
  void DeleteConfiguration(
    Configuration configuration) {
    std::vector<int> *heads = static_cast<std::vector<int>*>(configuration);
    delete heads;
  }

  // Create configuration.
  Configuration CreateConfiguration() {
    std::vector<int>* heads = new std::vector<int>(length_);
    return static_cast<Configuration>(heads);
  }

 public:
  void Initialize(int length, const std::vector<Arc*> &arcs) {
    length_ = length;
    index_arcs_.assign(length, std::vector<int>(length, -1));
    for (int k = 0; k < arcs.size(); ++k) {
      int h = arcs[k]->head();
      int m = arcs[k]->modifier();
      index_arcs_[h][m] = k;
    }
  }

 private:
  void RunChuLiuEdmondsIteration(std::vector<bool> *disabled,
                                 std::vector<std::vector<int> > *candidate_heads,
                                 std::vector<std::vector<double> >
                                 *candidate_scores,
                                 std::vector<int> *heads,
                                 double *value);
 protected:
  int length_; // Sentence length (including root symbol).
  std::vector<std::vector<int> > index_arcs_;
};

} // namespace AD3

#endif // FACTOR_TREE_H_
