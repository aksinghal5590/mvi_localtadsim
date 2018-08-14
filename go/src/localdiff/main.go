package main

import (
	"hicutil"
	"flag"
	"strings"
	"fmt"
	"os"
	"bufio"
	"strconv"
	"sort"
	"math"
	"math/rand"
	"runtime"
	"runtime/pprof"
	"sync"
	"time"
	)

type bdyvi struct {
	start int
	end int
	vi float64
	pval float64
}

func sortAndPrint(bdyvis []bdyvi, isPrint bool) {
	sort.Slice(bdyvis, func(i,j int) bool {
		if bdyvis[i].start == bdyvis[j].start {return bdyvis[i].end < bdyvis[j].end} else {return bdyvis[i].start < bdyvis[j].start}
	})
	if isPrint {
		for i := 0; i < len(bdyvis); i++ {
			fmt.Printf("%d\t%d\t%6f\t%6f\n", bdyvis[i].start, bdyvis[i].end, bdyvis[i].vi, bdyvis[i].pval)
		}
	}
}


func main() {

	runtime.GOMAXPROCS(1)
	numCPU := runtime.NumCPU()
	// inputs should be 2 TAD files, a resolution parameter, and the output file name
	tadin := flag.String("tad", "", "comma-separated list of two TAD filenames or file patterns if optimizing gamma")
	gammain := flag.String("gamma", "", "if optimizing gamma, use 'opt,n' where n is the median TAD size to optimize for")
	res := flag.Int("res", 1, "resolution of Hi-C data")
	outfile := flag.String("o", "", "output filename")
	pcount := flag.Int("p", numCPU, "no of processes")
	cpuprof := flag.String("cpu", "", "write cpu profile to file")

	flag.Parse()

	if *cpuprof != "" {
		f,_ := os.Create(*cpuprof)
			pprof.StartCPUProfile(f)
			defer pprof.StopCPUProfile()
	}

	tadfilelist := strings.Split(*tadin, ",")

	var gammaopt bool
	medtadlen := 0.0
	var err error
	if len(*gammain) > 0 {
		gammaopt = true
		gammadata := strings.Split(*gammain,",")
		medtadlen,err = strconv.ParseFloat(gammadata[1], 64)
		if err != nil {
			fmt.Println("Error: couldn't convert median TAD length value to float, make sure to input i.e. '-gamma=opt,100' ")
			os.Exit(1)
		}
	} else {
		gammaopt = false
	}

	// read TAD files and process TAD lists to fill in non-TAD regions
	then := time.Now()
	tadlists := processTADLists(tadfilelist, res, gammaopt, medtadlen)
	duration := time.Since(then)
	fmt.Printf("processTADLists duration: %s\n", duration)
	// fmt.Println("done processing TAD lists, choosing optimal gamma")

	// calculate VI values at boundaries (using DP)
	then = time.Now()
	bdyvis := calcMultiVI(tadlists)
	sortAndPrint(bdyvis, true)
	//bdyvis = calcMultiVIPairwise(tadlists)
	//sortAndPrint(bdyvis, true)
	duration = time.Since(then)
	fmt.Printf("calcVIatBdys duration: %s\n", duration)

	// calculate all p-values, select significant points
	convergencecondition := 1e-5
	numCPU = *pcount
	runtime.GOMAXPROCS(numCPU)
	then = time.Now()
	sigpts := calcAllPvals(tadlists, bdyvis, numCPU, convergencecondition)
	duration = time.Since(then)
	fmt.Printf("calcAllPvals duration: %s\n", duration)
	// fmt.Println("done calculating all p-values")

	// identify dominating points from significant ones
	runtime.GOMAXPROCS(1)
	then = time.Now()
	dompts := findDomPts(sigpts)
	duration = time.Since(then)
	fmt.Printf("findDomPts duration: %s\n", duration)
	// fmt.Println("done finding dominating points")

	// save results to a file
	writeOutputToFile(dompts,outfile)
}


func processTADLists(tadfilelist []string, res *int, gammaopt bool, medtadlen float64) ([][][]int) {

	N := len(tadfilelist)

	tadlists := make([][][]int, N)
	chrlength := 0
	var gamma float64
	for i:=0; i < N; i++ {
		if gammaopt == true {
			tadlists[i],gamma = hicutil.ChooseGamma(medtadlen, tadfilelist[i], *res)
			_ = gamma
		} else {
			tadlists[i] = hicutil.ReadTADFile(tadfilelist[i], *res)
		}
		n := tadlists[i][len(tadlists[i])-1][1]
		if chrlength < n+1 {
			chrlength = n+1
		}
	}
	for i := 0; i < N; i++ {
		tadlists[i] = hicutil.FillinTADList(tadlists[i], chrlength)
	}

	return tadlists

}

func calcVIatBdysNaive(tadlists [][][]int) ([]bdyvi) {
	var bdyvilist []bdyvi
	var newbdyvi bdyvi
	for _,tadlist := range tadlists {
		for i,tadstart := range tadlist {
			for _,tadend := range tadlist[i:] {

				if hicutil.ContainsHangingTAD(tadlists[0], tadstart[0], tadend[1]) || hicutil.ContainsHangingTAD(tadlists[1], tadstart[0], tadend[1]) {
					continue
				}

				newbdyvi.start = tadstart[0]
				newbdyvi.end = tadend[1]
				n := tadend[1] - tadstart[0] + 1
				intvl1 := hicutil.ProcessIntervals(tadlists[0],tadstart[0],tadend[1])
				intvl2 := hicutil.ProcessIntervals(tadlists[1],tadstart[0],tadend[1])
				overlaps := hicutil.CalcOverlaps(intvl1,intvl2)
				clus1sizes := make([]int, len(intvl1))
				for c,clus := range intvl1 {
					clus1sizes[c] = clus[1]-clus[0]+1 }
				clus2sizes := make([]int, len(intvl2))
				for c,clus := range intvl2 {
					clus2sizes[c] = clus[1]-clus[0]+1 }
				condh1 := hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
				condh2 := hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
				newbdyvi.vi = (condh1 + condh2)/math.Log(float64(n))
				bdyvilist = append(bdyvilist, newbdyvi)
			}
		}
	}
	return bdyvilist
}

func transpose(a [][]int) [][]int {
	n := len(a)
	m := len(a[0])
	b := make([][]int, m)
	for j := 0; j < m; j++ {
		b[j] = make([]int, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			b[j][i] = a[i][j]
		}
	}
	return b
}

func calcMultiVIPairwise(tadlists [][][]int) []bdyvi {
	var bdyvilist []bdyvi
	var newbdyviList []bdyvi
	var resultviList []bdyvi

	for i, tadlist1 := range tadlists {
		for j, tadlist2 := range tadlists {
			if j <= i {
				continue
			}
			tempTadList := [][][]int{tadlist1, tadlist2}
			bdyList := calcVIatBdysNaive(tempTadList)
			for _, newbdyvi := range bdyList {
				bdyvilist = append(bdyvilist, newbdyvi)
			}
		}
	}
	N := len(tadlists)
	for i, bvi1 := range bdyvilist {
		copyExists := 0
		for j, bvi2 := range bdyvilist {
			if i == j {
				continue
			}
			if bvi2.start == bvi1.start && bvi2.end == bvi1.end {
				copyExists++
				bvi1.vi += bvi2.vi
			}
		}
		if copyExists == N - 1 {
			newbdyviList = append(newbdyviList, bvi1)
			newbdyviList[len(newbdyviList) - 1].vi = newbdyviList[len(newbdyviList) - 1].vi / (float64)(N)
		}
	}
	// remove duplicates
	sortAndPrint(newbdyviList, false)
	for i := 0; i < len(newbdyviList); {
		bvi1 := newbdyviList[i]
		for j := i + 1;; {
			if newbdyviList[j].start == bvi1.start && newbdyviList[j].end == bvi1.end {
				j++
				if j >= len(newbdyviList) {
					i = j
					break
				}
				continue
			}
			i = j
			break
		}
		resultviList = append(resultviList, bvi1)
	}
	return resultviList
}

func calcMultiVI(tadlists [][][]int) []bdyvi {
	var bdyvilist []bdyvi
	var newbdyviList []bdyvi
	var resultviList []bdyvi
	N := len(tadlists)
	overlapSet := hicutil.GetOverlapSet(tadlists)
	for _, tadlist := range tadlists {
		bdyList := calcVI(overlapSet, tadlist)
		for _, newbdyvi := range bdyList {
			bdyvilist = append(bdyvilist, newbdyvi)
		}
	}
	for i, bvi1 := range bdyvilist {
		copyCount := 1
		for j, bvi2 := range bdyvilist {
			if j == i {
				continue
			}
			if bvi2.start == bvi1.start && bvi2.end == bvi1.end {
				copyCount++
				bvi1.vi += bvi2.vi
			}
		}
		if copyCount == N {
			newbdyviList = append(newbdyviList, bvi1)
			newbdyviList[len(newbdyviList) - 1].vi = newbdyviList[len(newbdyviList) - 1].vi / (float64)(copyCount)
		}
	}
	// remove duplicates
	sortAndPrint(newbdyviList, false)
	for i := 0; i < len(newbdyviList); {
		bvi1 := newbdyviList[i]
		for j := i + 1;; {
			if j >= len(newbdyviList) {
				i = j
				break
			}
			if newbdyviList[j].start == bvi1.start && newbdyviList[j].end == bvi1.end {
				j++
				continue
			}
			i = j
			break
		}
		resultviList = append(resultviList, bvi1)
	}
	return resultviList
}

func calcVI(overlapSet [][]int, tadlist [][]int) []bdyvi {
	var bdyvilist []bdyvi
	var s_pos int
	var newbdyvi bdyvi
	for i, tadstart := range overlapSet {
		for j := 0; j < len(tadlist); j++ {
			if tadstart[0] == tadlist[j][0] {
				s_pos = j
				break
			}
		}
		for t, tadend := range overlapSet[i:] {
			if hicutil.ContainsHangingTAD(tadlist, tadstart[0], tadend[1]) {
				continue
			}
			n := tadend[1] - tadstart[0] + 1
			vi := 0.0
			j := s_pos
			for pos := i; pos <= i + t; pos++ {
				st := tadlist[j][0]
				end := tadlist[j][1]
				P_i := (float64)(end - st + 1)
				P_intersect := (float64)(overlapSet[pos][1] - overlapSet[pos][0] + 1)
				vi += (P_intersect * math.Log(P_i/P_intersect))/(float64)(n)
				if tadlist[j][1] == overlapSet[pos][1] {
					j++
				}
			}
			vi = vi/math.Log((float64)(n))
			newbdyvi.start = tadstart[0]
			newbdyvi.end = tadend[1]
			newbdyvi.vi = vi
			bdyvilist = append(bdyvilist, newbdyvi)
		}
	}
	return bdyvilist
}

var wg sync.WaitGroup
func worker(tadlists [][][]int, job []bdyvi, result *[]bdyvi, convcond float64) {
	defer wg.Done()
	concurrentRand := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i, querypt := range job {
		(*result)[i] = appendPval(tadlists, querypt, convcond, concurrentRand)
	}
}

func calcAllPvals(tadlists [][][]int, bdyvis []bdyvi, numCPU int, convcond float64) []bdyvi {

	var sigpts []bdyvi
	if len(bdyvis) < numCPU {
		numCPU = len(bdyvis)
	}
	bdyvis_pval := make([]bdyvi, len(bdyvis))
	allpvals := make([]float64, len(bdyvis))
	jobs := make([][]bdyvi, numCPU)
	results := make([][]bdyvi, numCPU)

	k := 0
	for i := 0; i < len(bdyvis); i++ {
		jobs[k] = append(jobs[k], bdyvis[i])
		k = (k + 1) % numCPU
	}

	for w, job := range jobs {
		wg.Add(1)
		results[w] = make([]bdyvi, len(job))
		go worker(tadlists, job, &results[w], convcond)
	}
	wg.Wait()

	i := 0
	for w := 0; w < numCPU; w++ {
		for _, res := range results[w] {
			bdyvis_pval[i] = res
			allpvals[i] = bdyvis_pval[i].pval
			i++
		}
	}

	bhidx := hicutil.MultHypTestBH(allpvals)
	if bhidx > -1 {
		sort.Slice(bdyvis_pval, func(i,j int) bool {return bdyvis_pval[i].pval < bdyvis_pval[j].pval})
		sigpts = bdyvis_pval[:bhidx+1]
	}
	return sigpts
}

func appendPval( tadlists [][][]int, querypt bdyvi, convcond float64, r *rand.Rand) (bdyvi) {

	N := len(tadlists)
	intvls := make([][][]int, N)
	for i := 0; i < N; i++ {
		intvls[i] = hicutil.ProcessIntervals(tadlists[i], querypt.start, querypt.end)
	}
	n := querypt.end - querypt.start + 1
	p := hicutil.CalcPval(intvls, n, querypt.vi, convcond, r)
	isSize1 := false
	for i := 0; i < N; i++ {
		if len(intvls[i]) == 1 {
			isSize1 = true
			break
		}
	}
	if p < 0.05 && isSize1 == true {
		fmt.Println(p)
		fmt.Println(n,querypt.vi)
		fmt.Println(querypt.start, querypt.end)
		os.Exit(1)
	}
	querypt.pval = p
	return querypt
}

func findDomPts(sigpts []bdyvi) []bdyvi {

	var dompts []bdyvi
	for i,querypt := range sigpts {
		isdom := true
		for j,comppt := range sigpts {
			if i == j { continue }
			if comppt.start == querypt.start && comppt.end > querypt.end && comppt.vi < querypt.vi {
				isdom = false
				break
			}
			if comppt.start < querypt.start && comppt.end == querypt.end && comppt.vi < querypt.vi {
				isdom = false
				break
			}
			if comppt.start <= querypt.start { continue }
			if comppt.end >= querypt.end { continue }
			if comppt.vi < querypt.vi {
				isdom = false
				break
			}
		}
		if isdom == true {
			dompts = append(dompts, querypt)
		}
	}
	//if there are multiple dominating, significant points that start or end at the same place, remove the smaller interval
	var toremove []int
	for i,dompt1 := range dompts {
		for j,dompt2 := range dompts {
			if i == j {continue}
			if dompt1.start == dompt2.start && dompt1.end < dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start > dompt2.start && dompt1.end == dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start == dompt2.start && dompt1.end == dompt2.end && i < j {
				toremove = append(toremove, i)
				break
			}
		}
	}
	sort.Ints(toremove)
	for i,j := range toremove {
		dompts = append(dompts[:j-i], dompts[j-i+1:]...)
	}
	sort.Slice(dompts, func(i,j int) bool {return dompts[i].start < dompts[j].start})
	return dompts
}

func writeOutputToFile(domsigpts []bdyvi, outfile *string) {

    //write values to file
    f,err := os.Create(*outfile)
    if err != nil {
            panic(err)
    }
    //defer f.Close()

    w := bufio.NewWriter(f)
    labelline := []string{"start", "end", "VI", "p-value"}
    //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
    //fmt.Println(line1)
    fmt.Fprintf(w,line1+"\n")

    for _,vals := range domsigpts {
            strvals := make([]string, 4)
	strvals[0] = strconv.Itoa(vals.start)
	strvals[1] = strconv.Itoa(vals.end)
	strvals[2] = strconv.FormatFloat(vals.vi,'g',6,64)
	strvals[3] = strconv.FormatFloat(vals.pval,'g',6,64)
            newline := strings.Join(strvals, "\t")
            fmt.Fprintf(w,newline+"\n")
    }
    w.Flush()
    f.Close()
	//fmt.Println("Wrote output values to", *outfile)
}

