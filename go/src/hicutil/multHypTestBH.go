package hicutil

import (
	"sort"
	//"fmt"
)


func MultHypTestBH(pvals []float64) int {
	// multiple hypothesis correction through Benjamini-Hochberg procedure, at level 0.05

	alpha := 0.05
	numtests := len(pvals)
	// order lists by p-values
	sort.Float64s(pvals)
	// define thresholds
	thresh := make([]float64, numtests)
	for i := 0; i < numtests; i++ {
		thresh[i] = float64(i)*alpha/float64(numtests)
	}
	// find greatest index where pval < thresh
	imax := -1
	for i := numtests-1; i >= 0; i-- {
		if pvals[i] < thresh[i] {
			imax = i
			break
		}
	}
	return imax
}
