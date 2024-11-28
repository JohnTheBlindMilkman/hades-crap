/**
 * @file TrackSmearer.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Class implementing the momentum smearing of TrackCandidate objects. A word of caution: I do not change ALL of the TracCandidate variables, just Px,Py,Pz, and E (for femtoscopy)
 * @version 0.1.0
 * @date 2024-10-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TrackSmearer_hxx
    #define TrackSmearer_hxx

    #include <memory>
    #include <random>

    #include "TF1.h"

    #include "TrackCandidate.hxx"

    namespace Correction
    {
        class TrackSmearer
        {
            public:
                /**
                 * @brief Enum defining the available function types
                 * 
                 */
                enum class FunctionType{MomentumMean,MomentumStdev,PhiMean,PhiStdev,ThetaMean,ThetaStdev};

                /**
                 * @brief Construct a new Track Smearer object
                 * 
                 */
                TrackSmearer() = delete;
                /**
                 * @brief Construct a new Track Smearer object
                 * 
                 * @param momMu function for mean parameter of momentum
                 * @param momSigma function for standard deviation parameter of momentum
                 * @param phiMu function for mean parameter of azimuthal angle
                 * @param phiSigma function for standard deviation parameter of azimuthal angle
                 * @param thetaMu function for mean parameter of polar angle
                 * @param thetaSigma function for standard deviation parameter of polar angle
                 */
                TrackSmearer(std::unique_ptr<TF1> &&momMu,std::unique_ptr<TF1> &&momSigma,std::unique_ptr<TF1> &&phiMu,std::unique_ptr<TF1> &&phiSigma,std::unique_ptr<TF1> &&thetaMu,std::unique_ptr<TF1> &&thetaSigma);
                /**
                 * @brief Destroy the Track Smearer object
                 * 
                 */
                ~TrackSmearer() = default;
                /**
                 * @brief Construct a new Track Smearer object
                 * 
                 */
                TrackSmearer(const TrackSmearer&) = delete;
                /**
                 * @brief Copy-assign a new Track Smearer object
                 * 
                 * @return TrackSmearer& 
                 */
                TrackSmearer& operator=(const TrackSmearer&) = delete;
                /**
                 * @brief Construct a new Track Smearer object
                 * 
                 */
                TrackSmearer(TrackSmearer&&) = delete;
                /**
                 * @brief Move-assign a new Track Smearer object
                 * 
                 * @return TrackSmearer& 
                 */
                TrackSmearer& operator=(TrackSmearer&&) = delete;
                /**
                 * @brief Set the Fit Function object
                 * 
                 * @param func 
                 * @param type 
                 */
                void SetFitFunction(std::unique_ptr<TF1> &&func, FunctionType type) noexcept;
                /**
                 * @brief Function for smearing the momenta of given TrackCandidate object
                 * 
                 * @param cand 
                 */
                void SmearMomenta(Selection::TrackCandidate &cand);

            private:
                std::unique_ptr<TF1> m_meanMomentumFit, m_meanPhiFit, m_meanThetaFit, m_stdevMomentumFit, m_stdevPhiFit, m_stdevThetaFit;
                std::random_device m_randDevice;
                std::mt19937 m_messerneTwisterEngine;
                std::normal_distribution<double> m_momentumDistribution, m_phiDistribution, m_thetaDistribution;
                float m_trueMomentum, m_truePhi, m_trueTheta, m_smearedMomentum, m_smearedPhi, m_smearedTheta, m_momentumX, m_momentumY, m_momentumZ, m_energy;

        };
        
        TrackSmearer::TrackSmearer(std::unique_ptr<TF1> &&momMu,std::unique_ptr<TF1> &&momSigma,std::unique_ptr<TF1> &&phiMu,std::unique_ptr<TF1> &&phiSigma,std::unique_ptr<TF1> &&thetaMu,std::unique_ptr<TF1> &&thetaSigma) :
        m_meanMomentumFit(std::move(momMu)), m_meanPhiFit(std::move(phiMu)), m_meanThetaFit(std::move(thetaMu)), 
        m_stdevMomentumFit(std::move(momSigma)), m_stdevPhiFit(std::move(phiSigma)), m_stdevThetaFit(std::move(thetaSigma)),
        m_randDevice(), m_messerneTwisterEngine(m_randDevice()), m_momentumDistribution(0,1), m_phiDistribution(0,1), m_thetaDistribution(0,1),
        m_trueMomentum(0.), m_truePhi(0.), m_trueTheta(0.), m_smearedMomentum(0.), m_smearedPhi(0.), m_smearedTheta(0.), m_momentumX(0.),
        m_momentumY(0.), m_momentumZ(0.), m_energy(0.)
        {
        }

        void TrackSmearer::SetFitFunction(std::unique_ptr<TF1> &&func, FunctionType type) noexcept
        {
            switch (type)
            {
                case FunctionType::MomentumMean:
                    m_meanMomentumFit = std::move(func);
                    break;
                case FunctionType::MomentumStdev:
                    m_stdevMomentumFit = std::move(func);
                    break;
                case FunctionType::PhiMean:
                    m_meanPhiFit = std::move(func);
                    break;
                case FunctionType::PhiStdev:
                    m_stdevPhiFit = std::move(func);
                    break;
                case FunctionType::ThetaMean:
                    m_meanThetaFit = std::move(func);
                    break;
                case FunctionType::ThetaStdev:
                    m_stdevThetaFit = std::move(func);
                    break;

                
                default:
                    break;
            }
        }

        void TrackSmearer::SmearMomenta(Selection::TrackCandidate &cand)
        {
            m_trueMomentum = cand.GetP();
            m_truePhi = cand.GetPhi();
            m_trueTheta = cand.GetTheta();

            m_momentumDistribution = std::normal_distribution<double>(m_meanMomentumFit->Eval(m_trueMomentum),m_stdevMomentumFit->Eval(m_trueMomentum));
            m_phiDistribution = std::normal_distribution<double>(m_meanPhiFit->Eval(m_truePhi),m_stdevPhiFit->Eval(m_truePhi));
            m_thetaDistribution = std::normal_distribution<double>(m_meanThetaFit->Eval(m_trueTheta),m_stdevThetaFit->Eval(m_trueTheta));

            m_smearedMomentum = m_trueMomentum + m_momentumDistribution(m_messerneTwisterEngine);
            m_smearedPhi = m_truePhi + m_phiDistribution(m_messerneTwisterEngine);
            m_smearedTheta = m_trueTheta + m_thetaDistribution(m_messerneTwisterEngine);
            
            m_momentumX = m_smearedMomentum * std::cos(m_smearedPhi) * std::sin(m_smearedTheta);
            m_momentumY = m_smearedMomentum * std::sin(m_smearedPhi) * std::sin(m_smearedTheta);
            m_momentumZ = m_smearedMomentum * std::cos(m_smearedTheta);
            // E^2 = (pc)^2 + (m_0c^2)^2
            m_energy = std::sqrt(m_momentumX * m_momentumX + m_momentumY * m_momentumY + m_momentumZ * m_momentumZ + cand.GetM2());

            cand.SetMomentum(m_momentumX,m_momentumY,m_momentumZ,m_energy);
        }
    }

#endif