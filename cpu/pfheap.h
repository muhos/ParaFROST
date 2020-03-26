/***********************************************************************[pfheap.h]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Technische Universiteit Eindhoven (TU/e).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

#ifndef __HEAP_
#define __HEAP_

#include "pfdefs.h"

typedef struct {
	double var_inc, var_decay;
}VAR_PARAM;
class VAR_HEAP {
	uVector1D heap;
	vector1D indices;  // var index in heap
	double* activity;
	uint32 alloc_sz;
	VAR_PARAM var_param;
	static inline int left(int i) { return (i << 1) + 1; }
	static inline int right(int i) { return (i + 1) << 1; }
	static inline int parent(int i) { return (i - 1) >> 1; }

	void bubbleUp(int i)
	{
		uint32   x = heap[i];
		int p = parent(i);
		while (i != 0 && activity[x] > activity[heap[p]]) {
			heap[i] = heap[p];
			indices[heap[p]] = i;
			i = p;
			p = parent(p);
		}
		heap[i] = x;
		indices[x] = i;
	}

	void bubbleDown(int i)
	{
		uint32 x = heap[i];
		while (left(i) < heap.size()) {
			int child = (right(i) < heap.size() && activity[heap[right(i)]] > activity[heap[left(i)]]) ? right(i) : left(i);
			if (activity[heap[child]] <= activity[x]) break;
			heap[i] = heap[child];
			indices[heap[i]] = i;
			i = child;
		}
		heap[i] = x;
		indices[x] = i;
	}

public:
	VAR_HEAP() {
		var_param = { 0, 0 };
		activity = NULL;
		alloc_sz = 0;
	}
	~VAR_HEAP() {
		delete[] activity;
		activity = NULL;
		heap.clear(true);
		indices.clear(true);
	}
	void allocMem(const uint32& sz) {
		alloc_sz = sz;
		assert(alloc_sz > 0);
		activity = new double[alloc_sz];
		indices.resize(alloc_sz);
	}

	void init(double var_inc, double var_decay) {
		assert(indices.size() <= (int)alloc_sz);
		heap.clear();
		for (uint32 i = 0; i < alloc_sz; i++) {
			activity[i] = 0.0;
			indices[i] = UNDEFINED;
		}
		var_param.var_inc = var_inc;
		var_param.var_decay = var_decay;
	}

	void build(uint32 nVars) {
		for (uint32 v = 0; v < nVars; v++) insert(v);
	}

	void incVarDecay(double rate) {
		var_param.var_decay += rate;
	}

	double varActivity(const uint32& v) {
		return *(activity + v);
	}

	void varBumpAct(const uint32& v) {
		if ((activity[v] += var_param.var_inc) > 1e100) {
			for (uint32 i = 0; i < nOrgVars(); i++) activity[i] *= 1e-100; // rescale activity
			var_param.var_inc *= 1e-100;
		}
		if (has(v)) bubbleUp(indices[v]);
	}

	void varBumpAct(const uint32& v, double norm_act) {
		if ((activity[v] += norm_act) > 1e100) {
			for (uint32 i = 0; i < nOrgVars(); i++) activity[i] *= 1e-100; // rescale activity
			var_param.var_inc *= 1e-100;
		}
		if (has(v)) bubbleUp(indices[v]);
	}

	void VarDecayAct(void) {
		var_param.var_inc *= (1.0 / var_param.var_decay);
	}

	double VarDecay(void) {
		return var_param.var_decay;
	}

	int size() {
		return (int)heap.size();
	}

	bool empty() {
		return heap.size() == 0;
	}

	bool has(const uint32& v) {
		return indices[v] >= 0;
	}

	void insert(const uint32& v)
	{
		assert(v < (uint32)indices.size());
		assert(!has(v));
		indices[v] = (int)heap.size();
		heap.push(v);
		bubbleUp(indices[v]);
	}

	void remove(const uint32& v)
	{
		assert(has(v));
		int v_pos = indices[v];
		indices[v] = UNDEFINED;
		int last_idx = heap.size() - 1;
		if (v_pos < last_idx) {
			heap[v_pos] = heap[last_idx];
			indices[heap[v_pos]] = v_pos;
			heap.pop();
			bubbleDown(v_pos);
		}
		else heap.pop();
	}

	uint32 removeMin()
	{
		uint32 x = heap[0];
		heap[0] = heap[heap.size() - 1];
		indices[heap[0]] = 0;
		indices[x] = UNDEFINED;
		heap.pop();
		if (heap.size() > 1) bubbleDown(0);
		return x;
	}

	double maxActivity() {
		return activity[heap[0]] == 0.0 ? var_param.var_inc : activity[heap[0]];
	}

	double varInc() {
		return var_param.var_inc;
	}

	uint32* h_ptr() {
		return heap;
	}

	double* act_ptr() {
		return activity;
	}

	void clear() {
		for (int i = 0; i < heap.size(); i++) indices[heap[i]] = UNDEFINED;
		heap.clear(true);
		indices.clear(true);
	}

	void print() { // print heap content with associated activity
		for (int i = 0; i < heap.size(); i++) {
			printf(" heap[%d] = %d, a = %f\n", i, heap[i], activity[heap[i]]);
		}
	}
};

#endif