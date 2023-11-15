#ifndef UTILS_H
#define UTILS_H

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TColor.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <type_traits>

namespace util {

enum LogLevel {
    pERROR,
    pWARNING,
    pINFO
};

class Logger {
public:
    Logger(const std::string& tag, LogLevel level) : tag_(tag), level_(level) {}

    template <typename T>
    Logger& operator<<(const T& msg) {
        std::ostringstream oss;
        oss << msg;
        message_ += oss.str();
        return *this;
    }

    ~Logger() {
        std::string levelStr;
        switch (level_) {
            case pERROR: levelStr = "ERROR"; break;
            case pWARNING: levelStr = "WARNING"; break;
            case pINFO: levelStr = "INFO"; break;
        }
        std::cout << "[" << tag_ << "] " << levelStr << ": " << message_ << std::endl;
    }

private:
    std::string tag_;
    LogLevel level_;
    std::string message_;
};

#define LOG(tag, level) Logger(tag, level)

void Die(int exitCode) {
    std::exit(exitCode);
}

// Start from Steve Dennis
TMatrixD CovarianceMatrixToCorrelationMatrix(const TMatrixD & cov_m) {
    int n = cov_m.GetNrows();
    TMatrixD corr_m(n,n);
    for (int ii = 0 ; ii < n; ii++) {
        for (int jj = 0 ; jj < n; jj++) {
            corr_m[ii][jj] = cov_m[ii][jj] / sqrt(cov_m[ii][ii] * cov_m[jj][jj]);
        }
    }
    return corr_m;
}

TVectorD CovarianceMatrixToErrorVector(const TMatrixD & cov_m) {
    int n = cov_m.GetNrows();
    TVectorD vec(n);
    for (int ii = 0 ; ii < n; ii++) {
        vec[ii] = std::sqrt(cov_m[ii][ii]);
    }
    return vec;
}

TMatrixD CovarianceMatrixToFractionalCovarianceMatrix(const TMatrixD & cov_m, const TVectorD & mean) {
    int n = cov_m.GetNrows();
    TMatrixD frac_cov(n,n);
    for (int ii = 0 ; ii < n; ii++) {
        for (int jj = 0 ; jj < n; jj++) {
            frac_cov[ii][jj] = cov_m[ii][jj] / (mean[ii]*mean[jj]);
        }
    }
    return frac_cov;
}

TVectorD ErrorVectorToFractionalErrorVector(const TVectorD & total_err, const TVectorD & mean) {
    int n = total_err.GetNrows();
    TVectorD frac_err(n);
    for (int ii = 0 ; ii < n; ii++) {
        frac_err[ii] = total_err[ii] / mean[ii];
    }
    return frac_err;
}

TMatrixD CorrelationAndErrorToCovarianceMatrix(const TMatrixD & corr, const TVectorD & one_sigma) {
    int n = corr.GetNrows();
    if (n!=one_sigma.GetNrows()) {
        LOG("UtlNumStat",pERROR)<<"Matrix and vector dimensions do not match.";
        util::Die(-1);
        return TMatrixD();
    }
    TMatrixD cov_m(n,n);
    for (int ii = 0 ; ii < n; ii++) {
        for (int jj = 0 ; jj < n; jj++) {
            cov_m[ii][jj] = corr[ii][jj] * one_sigma[ii] * one_sigma[jj];
        }
    }
    return cov_m;
}

TMatrixD CountsToConfusionMatrix(const TMatrixD & event_counts, const char* normalization = "row") {
    TMatrixD confusion_mat(event_counts);
    if (strcmp(normalization, "row") == 0) {
        for (int i = 0; i < event_counts.GetNrows(); i++) {
            double row_sum = 0;
            for (int j = 0; j < event_counts.GetNcols(); j++) {
                row_sum += event_counts(i, j);
            }
            for (int j = 0; j < event_counts.GetNcols(); j++) {
                confusion_mat(i, j) /= row_sum;
            }
        }
    } else if (strcmp(normalization, "column") == 0) {
        for (int j = 0; j < event_counts.GetNcols(); j++) {
            double col_sum = 0;
            for (int i = 0; i < event_counts.GetNrows(); i++) {
                col_sum += event_counts(i, j);
            }
            for (int i = 0; i < event_counts.GetNrows(); i++) {
                confusion_mat(i, j) /= col_sum;
            }
        }
    } else
    {
        LOG("UtlNumStat",pERROR)<<"Normalization type not recognized.";
        util::Die(-1);
    }
    return confusion_mat;
}

int CreateRedToBlueColorPalette(int n_shades) {
    double reds  [3] = {1.,1.,0};
    double blues [3] = {0.,1.,1};
    double greens[3] = {0.,1.,0};
    
    double stops [3] = {0.0,0.5,1.0};

    int pal_number = TColor::CreateGradientColorTable(3,stops,reds,greens,blues,n_shades);
    return pal_number;
}

int CreateWhiteToBlueColorPalette(int n_shades) {
    double reds  [2] = {1.,0};
    double blues [2] = {1.,1};
    double greens[2] = {1.,0};
    
    double stops [2] = {0.0,1.0};

    int pal_number = TColor::CreateGradientColorTable(3,stops,reds,greens,blues,n_shades);
    return pal_number;
}

// End from Steve Dennis


TMatrixD TH2DToTMatrixD(const TH2D & hist) {
    int nX = hist.GetNbinsX();
    int nY = hist.GetNbinsY();
    TMatrixD mat(nX, nY);
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            mat[i][j] = hist.GetBinContent(i+1, j+1);
        }
    }
    return mat;
}

TH2D TMatrixDToTH2D(const TMatrixD & mat, const char* name, const char* title, double xlow, double xup, double ylow, double yup) {
    int nX = mat.GetNrows();
    int nY = mat.GetNcols();
    TH2D hist(name, title, nX, xlow, xup, nY, ylow, yup);
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            hist.SetBinContent(i+1, j+1, mat[i][j]);
        }
    }
    return hist;
}

template <typename T>
void PrintMatrix(const T & obj) {
    int nX = obj.GetNbinsX();
    int nY = obj.GetNbinsY() ? obj.GetNbinsY() : 1; // If GetNbinsY() is not available, set nY to 1

    // Check if T is a TMatrixD
    const int isTMatrixD = std::is_same<T, TMatrixD>::value;

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            double value = nY == 1 ? obj.GetBinContent(i+1) : obj.GetBinContent(i+1, j+1);
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace util

#endif // UTILS_H