package hicutil

import (
//	"fmt"
//	"os"
)

func ContainsHangingTAD(tadlist [][]int, start int, end int) bool {

	for _,tad := range tadlist {
		if start > tad[1] {
			continue
		}
		if start >= tad[0] && start <= tad[1] { // interval starts in current TAD
			if start >= (int)(1.0 * (float64)(tad[1] - tad[0])/2.0) + tad[0] {
				return true
			}
		}
		if end >= tad[0] && end <= tad[1] { // interval ends in current TAD
			if end <= (int)(1.0 * (float64)(tad[1] - tad[0])/2.0) + tad[0] {
				return true
			}
		} else if end <= tad[0] { // past end of interval
			break
		}

	}
	return false
}

/*func ContainsHangingTAD(tadlist [][]int, start int, end int) bool {

	for _,tad := range tadlist {
		if start > tad[1] {
			continue
		}
		if start >= tad[0] && start <= tad[1] { // interval starts in current TAD
			if start != tad[0] {
				return true
			}
		}
		if end >= tad[0] && end <= tad[1] { // interval ends in current TAD
			if end != tad[1] {
				return true
			}
		} else if end <= tad[0] { // past end of interval
			break
		}
	}
	return false
}*/
