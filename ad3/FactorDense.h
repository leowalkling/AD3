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

#ifndef FACTOR_DENSE_H_
#define FACTOR_DENSE_H_

#include "GenericFactor.h"
#include "MultiVariable.h"

namespace AD3 {

class FactorDense : public GenericFactor {
 public:
  FactorDense() {}
  virtual ~FactorDense() { ClearActiveSet(); }

  int type() { return FactorTypes::FACTOR_MULTI_DENSE; }

  // Print as a string.
  void Print(std::ostream& stream) {
    stream << "DENSE";
    Factor::Print(stream);

    // Write the number of multi-variables.
    stream << " " << multi_variables_.size();

    // Write the number of states for each multi-variable.
    for (int k = 0; k < multi_variables_.size(); ++k) {
      int num_states = multi_variables_[k]->GetNumStates();
      stream << " " << num_states;
    }

    // Write the additional log-potentials.
    int num_configurations = GetNumConfigurations();
    for (int index = 0; index < num_configurations; ++index) {
        stream << " " << std::setprecision(9) << additional_log_potentials_[index];
    }
      
    stream << std::endl;
  }

  // Compute the score of a given assignment.
  void Maximize(const std::vector<double> &variable_log_potentials,
                const std::vector<double> &additional_log_potentials,
                Configuration &configuration,
                double *value) {
    std::vector<int> *states = static_cast<std::vector<int>*>(configuration);
    int best = -1;
    *value = -1e12;
    for (int index = 0;
         index < additional_log_potentials.size();
         ++index) {
      double score = additional_log_potentials[index];
      GetConfigurationStates(index, states);
      for (int i = 0; i < states->size(); ++i) {
        int variable_index = GetVariableIndex(i, (*states)[i]);
        score += variable_log_potentials[variable_index];
      }
      if (best < 0 || score > *value) {
        best = index;
        *value = score;
      }
    }
    assert(best >= 0);
    GetConfigurationStates(best, states);
  }

  // Compute the score of a given assignment.
  void Evaluate(const std::vector<double> &variable_log_potentials,
                const std::vector<double> &additional_log_potentials,
                const Configuration configuration,
                double *value) {
    const std::vector<int>* states =
        static_cast<const std::vector<int>*>(configuration);
    *value = 0.0;
    int offset_states = 0;
    for (int i = 0; i < states->size(); ++i) {
      int state = (*states)[i];
      *value += variable_log_potentials[offset_states + state];
      offset_states = variable_offsets_[i]; 
    }
    int index = GetConfigurationIndex(*states);
    *value += additional_log_potentials[index];
  }

  // Given a configuration with a probability (weight), 
  // increment the vectors of variable and additional posteriors.
  void UpdateMarginalsFromConfiguration(
    const Configuration &configuration,
    double weight,
    std::vector<double> *variable_posteriors,
    std::vector<double> *additional_posteriors) {
    const std::vector<int> *states =
        static_cast<const std::vector<int>*>(configuration);
    int offset_states = 0;
    for (int i = 0; i < states->size(); ++i) {
      int state = (*states)[i];
      (*variable_posteriors)[offset_states + state] += weight;
      offset_states = variable_offsets_[i]; 
    }
    int index = GetConfigurationIndex(*states);
    (*additional_posteriors)[index] += weight;
  }
  
  // Count how many common values two configurations have.
  int CountCommonValues(const Configuration &configuration1,
                        const Configuration &configuration2) {
    const std::vector<int> *states1 =
        static_cast<const std::vector<int>*>(configuration1);
    const std::vector<int> *states2 =
        static_cast<const std::vector<int>*>(configuration2);
    assert(states1->size() == states2->size());
    int count = 0;
    for (int i = 0; i < states1->size(); ++i) {
      if ((*states1)[i] == (*states2)[i]) ++count;
    }
    return count;
  }

  // Check if two configurations are the same.
  bool SameConfiguration(
    const Configuration &configuration1,
    const Configuration &configuration2) {
    const std::vector<int> *states1 = static_cast<const std::vector<int>*>(configuration1);
    const std::vector<int> *states2 = static_cast<const std::vector<int>*>(configuration2);
    assert(states1->size() == states2->size());
    for (int i = 0; i < states1->size(); ++i) {
      if ((*states1)[i] != (*states2)[i]) return false;
    }
    return true;
  }

  // Delete configuration.
  void DeleteConfiguration(
    Configuration configuration) {
    std::vector<int> *states = static_cast<std::vector<int>*>(configuration);
    delete states;
  }

  Configuration CreateConfiguration() {
    int length = multi_variables_.size();
    std::vector<int>* states = new std::vector<int>(length, -1);
    return static_cast<Configuration>(states); 
  }

 public:
  // Initialize the factor and build internal structure. 
  // num_states contains the number of states for each multi-variable linked to 
  // the factor.
  // Note: the variables and the additional log-potentials must be ordered
  // properly.
  void Initialize(const std::vector<MultiVariable*> &multi_variables) {
    multi_variables_ = multi_variables;

    // Build offsets.
    variable_offsets_.resize(multi_variables_.size());
    variable_offsets_[0] = multi_variables_[0]->GetNumStates();
    for (int i = 1; i < multi_variables_.size(); ++i) {
      variable_offsets_[i] = variable_offsets_[i-1] +
        multi_variables_[i]->GetNumStates();
    }
  }  

  int GetNumMultiVariables() { return multi_variables_.size(); }
  MultiVariable *GetMultiVariable(int i) { return multi_variables_[i]; }
  
  int GetNumConfigurations() {
    int  num_configurations = 1;
    for (int i = 0; i < multi_variables_.size(); ++i) {
      num_configurations *= multi_variables_[i]->GetNumStates();
    }
    return num_configurations;    
  }
 
  // Find the configuration index given the array of states.
  int GetConfigurationIndex(const std::vector<int>& states) {
    int index = states[0];
    for (int i = 1; i < states.size(); ++i) {
      index *= multi_variables_[i]->GetNumStates();
      index += states[i];
    }
    return index;
  }

  // Find the array of states corresponding to a configuration index.
  void GetConfigurationStates(int index, std::vector<int> *states) {
    int tmp = 1;
    for (int i = 1; i < states->size(); ++i) {
      tmp *= multi_variables_[i]->GetNumStates();
    }
    (*states)[0] = index / tmp;
    for (int i = 1; i < states->size(); ++i) {
      index = index % tmp;
      tmp /= multi_variables_[i]->GetNumStates();
      (*states)[i] = index / tmp;
    }
  }

  int GetVariableIndex(int i, int state) {
    return (i == 0)? state : variable_offsets_[i-1] + state;
  }

 private:
  // Multi-variables linked to this factor.
  std::vector<MultiVariable*> multi_variables_;
  // Offsets of the variables in the pool of values.
  std::vector<int> variable_offsets_;
};

} // namespace AD3

#endif
