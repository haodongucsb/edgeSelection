library(testthat)
library(edgeSelection)

test_check("edgeSelection")
# test_joint_matrix<-matrix(1,6,6)
# test_joint_matrix[2,6]<-test_joint_matrix[6,2]<-test_joint_matrix[4,5]<-test_joint_matrix[5,4]<-test_joint_matrix[5,6]<-test_joint_matrix[6,5]<-0
# test_that("test the joint method on NO2 dataset",
#           {
#             expect_equal(edge.selection(data = NO2, family = "joint", nbasis = 100)>0,test_joint_matrix)
#           })
