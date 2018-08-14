package hicutil

import (
	"math"
//	"fmt"
)

func findMinOverlap(tadlists [][][]int, posList []int) (int, int) {
	N := len(tadlists)
	maxStart := 0
	minEnd := math.MaxInt64
	for i := 0; i < N; i++ {
		if maxStart < tadlists[i][posList[i]][0] {
			maxStart = tadlists[i][posList[i]][0]
		}
		if minEnd > tadlists[i][posList[i]][1] {
			minEnd = tadlists[i][posList[i]][1]
		}
	}
	return maxStart, minEnd
}

func GetOverlapSet(tadlists [][][]int) [][]int {
	N := len(tadlists)
	posList := make([]int, N)
	sizeList := make([]int, N)
	var overlapSet [][]int
	for i := 0; i < N; i++ {
		posList[i] = 0
		sizeList[i] = len(tadlists[i])
	}
	for {
		isComplete := true
		for i := 0; i < N; i++ {
			if posList[i] <= sizeList[i] - 1 {
				isComplete = false
			}
		}
		if isComplete == true {
			break
		}

		maxStart, minEnd := findMinOverlap(tadlists, posList)
		temp := []int{maxStart, minEnd}
		overlapSet = append(overlapSet, temp)
		for i := 0; i < N; i++ {
			if minEnd >= tadlists[i][posList[i]][1] {
				posList[i]++
			}
		}
	}
	return overlapSet
}