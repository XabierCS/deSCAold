require(testthat)
require(deSCA)
pathX<-'/Users/xabier/Desktop/ipsych/reviseThesis/paper2/packagesR/deSCA/tests/data/array1.txt.gz'
print('Testing deSCA package, may the force be with you')
test_that('test SlideWindowMean',{
  v <-  c(2, 4, 6, 8, 10)
  b <- c(4, 6, 8, 9, 10)
  expect_equal(SlideWindowMean(v,2),b)
})

test_that('Test Reading raw files',
          { 
          Raw_file<-pathX
          Sample2<-ReadSample(RawFile = Raw_file,skip=10)
          
          expect_equal(sum(colnames(Sample2) %in% c("Log.R.Ratio", "B.Allele.Freq", "Chr")) ,3)
          }
)

test_that('Test Compute values from raw file',
          { 
            t1<-computeValues(path = pathX ,skip=10)
            dim(t1)
            
            expect_equal(dim(t1)  ,c(1,73))
          }
)



test_that('Test Simulate males',
          { 
            
            sample1 <- simulMales(10000,0,0,0.05,0.05,0.3,0.3)
            #1
            freq<-round(as.numeric(table(sample1$kariotype)/10000),3)
            expect_lt(freq[1]  ,0.06)
            expect_lt(freq[3]  ,0.06)
            expect_gt(freq[2]  ,freq[1])
            #2
            expect_equal(dim(sample1)  ,c(10000,3))
            #3
            expect_identical(names(sample1)  ,c('mX','mY','kariotype'))
          }
)



test_that('Test Simulate females',
          { 
            sample1 <- simulFemales(10000,0,0,0.05,0.05,0.3,0.3)
            #1
            freq<-round(as.numeric(table(sample1$kariotype)/10000),3)
            expect_lt(freq[1]  ,0.06)
            expect_lt(freq[3]  ,0.06)
            expect_gt(freq[2]  ,freq[1])
            #2
            expect_equal(dim(sample1)  ,c(10000,3))
            #3
            expect_identical(names(sample1)  ,c('mX','mA','kariotype'))
          }
)


test_that('Test male clusters',
          { 
            sample1 <- simulMales(10000,0,0,0.05,0.05,0.3,0.3)
            dataFrame2<-sample1[,c(1,2)]
            f1 <- optimizeClusterMales(dataFrame2, minPtsRange = c(12,15),interClusProb = 0.5)
            expect_identical(dim(f1[[2]]),as.integer(c(10000,8)))
            expect_identical(names(f1[[1]])[1:2],  c("cluster","minPts"))
            expect_identical(length(f1),as.integer(4))
            }
)



test_that('Test Female clusters',
          { 
            sample1 <- simulFemales(10000,0,0,0.05,0.05,0.3,0.3)
            dataFrame2<-sample1[,c(1,2)]
            f1 <- optimizeClusterFemales(dataFrame2, minPtsRange = c(12,15),interClusProb = 0.5)
            expect_identical(dim(f1[[2]]),as.integer(c(10000,8)))
            expect_identical(names(f1[[1]])[1:2],  c("cluster","minPts"))
            expect_identical(length(f1),as.integer(4))
          }
)


