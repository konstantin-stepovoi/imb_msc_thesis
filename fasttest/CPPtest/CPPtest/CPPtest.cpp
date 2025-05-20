#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <thread>
#include <mutex>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <optional>


float z_criterium(float mu1, float mu2, float sig1, float sig2, int n) {
    return (mu1 - mu2) / std::sqrt(sig1 / n + sig2 / n);
}

