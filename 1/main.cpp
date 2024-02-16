#include <gtest/gtest.h>
#include"solver.cpp"

TEST(Matrix_test, test_matrix) {
    vector solution=solver(TestMatrix.txt)
    ASSERT_NEAR(122, 32*result[0]+12*result[1], 0.0001);
    ASSERT_NEAR(324, 23*result[0]+34*result[1]+54*result[2], 0.0001);
    ASSERT_NEAR(127, 31*result[1]+76*result[2]+13*result[3], 0.0001);
    ASSERT_NEAR(64, 15*result[2]+54*result[3], 0.0001);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}