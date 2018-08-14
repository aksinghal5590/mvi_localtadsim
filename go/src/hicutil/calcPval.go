package hicutil

import (
	"math"
	"math/rand"
)

func CalcPval(intvls [][][]int, n int, vival float64, convcond float64, r *rand.Rand) float64 {

	N := len(intvls)
	pval := make([]float64, N)

	for i := 0; i < N; i++ {

		shuffnum := 0
		count := 0
		prevpvaldiffs := []float64{1000.0, 1000.0, 1000.0, 1000.0}
		prevpval := 1000.0
		keepshuffling := true

		tadlists := make([][][]int, N)
		copy(tadlists, intvls)

		clusSizes := make([]int, len(intvls[i]))
		for c,clus := range intvls[i] {
			clusSizes[c] = clus[1] - clus[0] + 1
		}
		for keepshuffling {
			shuffnum++
			newlist, _ := shuffledoms(clusSizes, r)
			copy(tadlists[i], newlist)
			overlapSet := GetOverlapSet(tadlists)
			j := 0
			vi := 0.0
			for k := 0; k < len(overlapSet); k++ {
				P_i := (float64)(newlist[j][1] - newlist[j][0] + 1)
				P_intersect := (float64)(overlapSet[k][1] - overlapSet[k][0] + 1)
				vi += (P_intersect * math.Log(P_i/P_intersect))/(float64)(n)
				if newlist[j][1] == overlapSet[k][1] {
					j++
				}
			}
			vi = vi/math.Log((float64)(n))
			if vi - vival < 1e-10 {
				count++
			}
			pval[i] = float64(count + 1)/float64(shuffnum + 1)
			prevpvaldiffs = append(prevpvaldiffs[1:], []float64{math.Abs(prevpval-pval[i])}...)
			// if last 5 p-values are within convergence condition, we're done shuffling
			keepshuffling = false
			for _,pvaldiff := range prevpvaldiffs {
				if pvaldiff > convcond {
					keepshuffling = true
				}
			}
			prevpval = pval[i]
		}
	}
	pv := 0.0
	for i := 0; i < N; i++ {
		pv += pval[i]
	}
	pv = pv/(float64)(N)
	return pv
}

func shuffledoms(clussizes []int, r *rand.Rand) ([][]int, []int) {

	permsizes := make([]int, len(clussizes))
	perm := r.Perm(len(clussizes))
	for j,v := range perm {
		permsizes[v] = clussizes[j]
	}
	// turn shuffled lists of lengths back into TAD lists
	newlist := make([][]int, len(clussizes))
	newlist[0] = []int{0,permsizes[0]-1}
	for j := 1; j < len(permsizes); j++ {
		tadstart := newlist[j-1][1]+1
		newlist[j] = []int{tadstart, tadstart+permsizes[j]-1}
	}

	return newlist, permsizes
}
