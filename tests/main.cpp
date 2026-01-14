#include <gtest/gtest.h>

// Najprostszy test — zawsze przejdzie
TEST(SampleTest, AlwaysPasses)
{
    EXPECT_EQ(2 + 2, 4);
}