/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan */

#ifndef COAL_BROADPHASE_DETAIL_INTERVALTREENODE_H
#define COAL_BROADPHASE_DETAIL_INTERVALTREENODE_H

#include "coal/broadphase/detail/simple_interval.h"
#include "coal/fwd.hh"

namespace coal {

namespace detail {

class COAL_DLLAPI IntervalTree;

/// @brief The node for interval tree
class COAL_DLLAPI IntervalTreeNode {
 public:
  friend class IntervalTree;

  /// @brief Create an empty node
  IntervalTreeNode();

  /// @brief Create an node storing the interval
  IntervalTreeNode(SimpleInterval* new_interval);

  ~IntervalTreeNode();

  /// @brief Print the interval node information: set left = invalid_node and
  /// right = root
  void print(IntervalTreeNode* left, IntervalTreeNode* right) const;

 protected:
  /// @brief interval stored in the node
  SimpleInterval* stored_interval;

  Scalar key;

  Scalar high;

  Scalar max_high;

  /// @brief red or black node: if red = false then the node is black
  bool red;

  IntervalTreeNode* left;

  IntervalTreeNode* right;

  IntervalTreeNode* parent;
};

}  // namespace detail
}  // namespace coal

#endif
