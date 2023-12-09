/*! \file    OCPFuncPVT.hpp
*   \brief   Functions for PVT in OCP
*   \author  Shizhe Li
*   \date    Jun/18/2023
*
*-----------------------------------------------------------------------------------
*  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
*  Released under the terms of the GNU Lesser General Public License 3.0 or later.  
*-----------------------------------------------------------------------------------
*/

#ifndef __OCPFUNCPVT_HEADER__ 
#define __OCPFUNCPVT_HEADER__  

// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"  
#include "UtilMath.hpp"
using namespace std;

/** @defgroup PVT ����
* @brief ���ڼ��������PVT����
* @details �������ڶ����͵�PVT���������Ǹ�����������PVT�������ܶȣ�ճ�ȵ�
* @{
*/

/** @brief ˮ�����ʱ����(PVTW)
* @details ����ˮ����ܶȣ�ճ�ȵ����ʣ������ǹ���ˮ��ѹ���ĺ���
*/
class OCP_PVTW : public OCPFuncTable
{
    public:
        //! Ĭ�Ϲ��캯��
        OCP_PVTW() = default; 

        //! ����PVT���
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoWin, const OCP_DBL& stdVwin);
        
        //! ����ˮ���Ħ��Ũ��
        OCP_DBL CalXiW(const OCP_DBL& P) const;  

        //! ����ˮ����ܶ�
        OCP_DBL CalRhoW(const OCP_DBL& P) const;

        //! ����ˮ��������ܶȣ�Ħ��Ũ�ȣ�ճ���Լ���صĵ�������
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

    protected:
        //! ����ˮ�������仯ϵ��
        OCP_DBL CalBw(const OCP_DBL& P) const;

        //! ����ˮ�������仯ϵ����ճ�ȵȵ���
        void CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;

    protected:
        //! ��̬��ˮ����ܶ�
        OCP_DBL stdRhoW;
        
        //! ��̬��ˮ���Ħ�����
        OCP_DBL stdVw;  
};

/** @brief ��ѹ���Ļ��͵�PVT���� 
* @details ����������ܶȣ�ճ�ȵ����ʣ������ǹ�������ѹ���ĺ���
*/
class OCP_PVCO : public OCPFuncTable  
{
    public:
        //! Ĭ�Ϲ��캯��
        OCP_PVCO() = default;

        //! ����PVCO���
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVoin, const OCP_DBL& stdVgin);
        
        //! ����������ܶ�
        OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! ���������Ħ��Ũ��
        OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! ���㱥��������ܶȣ�Ħ��Ũ�ȣ�ճ�ȣ����ͱȼ���Ӧ�ĵ���
        void CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                             OCP_DBL& rs, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP) const;
                             
        //! ���㲻����������ܶȣ�Ħ��Ũ�ȣ�ճ�ȼ���Ӧ�ĵ���
        void CalRhoXiMuDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, 
                           OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs) const;
                           
        //! ���㱥����������ͱ�
        OCP_DBL CalRs(const OCP_DBL& P) const;

    protected:
        //! ���㱥����������ͱȣ�����仯���ӣ�ճ�ȣ�����Ӧ�ĵ��� 
        void CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
                           OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP) const;
                           
        //! ���㲻�������������仯���ӣ�ճ�ȣ�����Ӧ�ĵ���
        void CalBoMuoDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu,  
                         OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs) const;
                         
    protected:
        //! ��̬�������ܶ�
        OCP_DBL stdRhoO;   
        
        //! ��̬�������ܶ�
        OCP_DBL stdRhoG;
        
        //! ��̬������Ħ�����
        OCP_DBL stdVo;
        
        //! ��̬������Ħ�����
        OCP_DBL stdVg;              
};

/** @brief ���������PVT����
* @details �������������ܶȣ�ճ�ȵ����ʣ������ǹ�������ѹ���ĺ���
*/
class OCP_PVDG : public OCPFuncTable
{
    public:
        //! Ĭ�Ϲ��캯��
        OCP_PVDG() = default;

        //! ����PVDG�����
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVgin);
                   
        //! ���������Ħ��Ũ��
        OCP_DBL CalXiG(const OCP_DBL& P) const;

        //! ����������ܶ�
        OCP_DBL CalRhoG(const OCP_DBL& P) const;

        //! ����������ܶȣ�Ħ��Ũ�ȣ�ճ�ȣ�����Ӧ�ĵ��� 
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                           
    protected:
        //! �������������仯���� 
        OCP_DBL CalBg(const OCP_DBL& P) const;

        //! �������������仯���ӣ�ճ�Ⱥ���Ӧ�ĵ���
        void CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, 
                         OCP_DBL& bP, OCP_DBL& muP) const;
                         
    protected:
        //! ��̬��������ܶ�
        OCP_DBL stdRhoG;
        
        //! ��̬�������Ħ�����
        OCP_DBL stdVg;
};

/** @brief ���͵�PVT����
* @details ������������ܶȣ�ճ�ȵ����ʣ������ǹ�������ѹ���ĺ���
*/  
class OCP_PVDO : public OCPFuncTable
{
    public:
        //! Ĭ�Ϲ��캯��
        OCP_PVDO() = default;

        //! ����PVDO�����
        virtual void Setup(const vector<vector<OCP_DBL>>& src, 
                           const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin);
                           
        //! ���������Ħ��Ũ��
        virtual OCP_DBL CalXiO(const OCP_DBL& P) const;

        //! ����������ܶ�
        virtual OCP_DBL CalRhoO(const OCP_DBL& P) const;

        //! ����������ܶȣ�Ħ��Ũ�ȣ�ճ�ȣ�����Ӧ�ĵ���
        virtual void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                                   OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                                   
    protected:
        //! �������������仯����  
        virtual OCP_DBL CalBo(const OCP_DBL& P) const;

        //! ���㲻�������������仯���ӣ�ճ�ȣ�����Ӧ�ĵ���
        virtual void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                                 OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                                 
    protected:
        //! ��̬��������ܶ�
        OCP_DBL stdRhoO;    
        
        //! ��̬�������Ħ�����
        OCP_DBL stdVo;     
};

/** @brief ���͵�PVT����(��ѹ��ϵ��Ϊ����)
* @details ������������ܶȣ�ճ�ȵ����ʣ������ǹ�������ѹ���ĺ�������ѹ��ϵ��Ϊ����
*/
class OCP_PVCDO : public OCP_PVDO
{
    public:
        //! Ĭ�Ϲ��캯��
        OCP_PVCDO() = default;

        //! ����PVCDO�����
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin) override;
                   
        //! ���������Ħ��Ũ��
        OCP_DBL CalXiO(const OCP_DBL& P) const override;

        //! ����������ܶ�
        OCP_DBL CalRhoO(const OCP_DBL& P) const override;

        //! ����������ܶȣ�Ħ��Ũ�ȣ�ճ�ȣ�����Ӧ�ĵ���
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const override;
                           
    protected:
        //! �������������仯����
        OCP_DBL CalBo(const OCP_DBL& P) const;

        //! ���㲻�������������仯���ӣ�ճ�ȣ�����Ӧ�ĵ��� 
        void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                         OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                         
    protected:
        //! �ο�ѹ��
        OCP_DBL Pref;
        
        //! �ο�ѹ���µ���������仯����
        OCP_DBL Bref;
        
        //! �ο�ѹ���µ������ѹ����ϵ��
        OCP_DBL Cb;
        
        //! �ο�ѹ���µ�����ճ��
        OCP_DBL muref;
        
        //! �ο�ѹ���µ�����ճ��ϵ��  
        OCP_DBL Cmu;
};

/** @brief һ��3άPVT��񣬱����������¶Ⱥ�ѹ���仯 
* @details ����������ܶȣ�ճ�ȵ����ʣ������ǹ�������ѹ�����¶ȵĺ���
*/
class OCP_PVT2 : public OCPFuncTable2
{
    public:
        //! ���������ܶ�
        OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T) const;

        //! �����ܶȣ�ճ�Ⱥ��ܽ��
        void CalRhoMuSol(const OCP_DBL& P, const OCP_DBL& T, 
                         OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol) const;
                         
        //! �����ܶȣ�ճ�Ⱥ��ܽ�Ⱥ���Ӧ�ĵ���
        void CalRhoMuSolDer(const OCP_DBL& P, const OCP_DBL& T, 
                            OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol,
                            OCP_DBL& rhoP, OCP_DBL& muP, OCP_DBL& solP) const;
};

// Typedefs for specialized PVT tables
typedef OCP_PVT2 OCP_PVTCO2; 
typedef OCP_PVT2 OCP_PVTH2O;

/** @brief Garciaw ģ�� 
* @details �����ܽ���CO2��ˮ���ܶ�
*/
class Garciaw
{
    public:
        //! ����Garciawģ�� 
        void Setup(const OCP_BOOL& flag);

        //! �Ƿ�ʹ�ø�ģ��
        auto IfUse() const;

        //! ����ˮ���ܶ�
        void CalRho(const OCP_DBL& T, const OCP_DBL& xGw, OCP_DBL& rhow) const;

        //! ����ˮ���ܶȺ���صĵ��� 
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP, 
                       OCP_DBL& rhow, OCP_DBL& rhowP, OCP_DBL& drhow_dxGw) const;
                       
        //! ����ˮ���ܶȺ���صĵ��� 
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP,  
                       const OCP_DBL& xGwT, OCP_DBL& rhow, OCP_DBL& rhowP, 
                       OCP_DBL& rhowT, OCP_DBL& drhow_dxGw) const;
                       
    protected:
        OCP_BOOL ifUse;      ///< �Ƿ�ʹ�ø�ģ��       
        const OCP_DBL MWCO2; ///< ������̼�������� 
};

/** @brief ճ�ȼ��������  
*/
class ViscosityParams 
{
    public:
        //! ���캯��
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin);

        //! ���캯�� 
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, 
                        const OCP_DBL* xin, const OCP_DBL* xiin);
                        
        //! ���캯��
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin,
                        const OCP_DBL* xin, const OCP_DBL* xiin,  
                        const OCP_DBL* xiPin, const OCP_DBL* xiTin,
                        const OCP_DBL* xixin);
                        
    public:
        const OCP_DBL* P;     ///< ѹ��
        const OCP_DBL* T;     ///< �¶�
        const OCP_DBL* x;     ///< ���Ħ������
        const OCP_DBL* xi;    ///< Ħ��Ũ��
        const OCP_DBL* xiP;   ///< Ħ��Ũ�ȶ�ѹ���ĵ���
        const OCP_DBL* xiT;   ///< Ħ��Ũ�ȶ��¶ȵĵ���
        const OCP_DBL* xix;   ///< Ħ��Ũ�ȶ�Ħ�������ĵ���
};

/** @brief ճ�ȼ��㷽�� 
*/
class ViscosityMethod 
{
    public:
        //! Ĭ�Ϲ��캯��
        ViscosityMethod() = default;

        //! ����ճ��
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp) = 0;

        //! ����ճ�Ⱥ���Ӧ�ĵ��� 
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp, 
                                     OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) = 0;
                                     
    protected:
        USI nc;                ///< �����
        vector<OCP_DBL> muc;   ///< ��ֵ�ճ��
        vector<OCP_DBL> mucP;  ///< ��ֵ�ճ�ȶ�ѹ���ĵ���
        vector<OCP_DBL> mucT;  ///< ��ֵ�ճ�ȶ��¶ȵĵ���
};

/** @brief ������ֵ�ѹ���¶�������ճ�ȱ����ճ�ȣ����Ի�Ϲ���
*/
class ViscosityMethod01 : public ViscosityMethod
{
    public:
        //! ���캯�� 
        ViscosityMethod01(const Table2& tab);

        //! ����ճ��
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! ����ճ�Ⱥ���Ӧ�ĵ���
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! ��ֵ�ѹ���¶�������ճ�ȱ�
        OCPTable2 viscTab;
};

/** @brief ������ֵ�ճ�ȹ�����������ճ�ȣ����Ի�Ϲ���
*/ 
class ViscosityMethod02 : public ViscosityMethod
{
    public:
        //! ���캯��
        ViscosityMethod02(const vector<OCP_DBL>& av, const vector<OCP_DBL>& bv);

        //! ����ճ�� 
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! ����ճ�Ⱥ���Ӧ�ĵ���
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! ��ֵ�ճ�ȹ������� 
        vector<OCP_DBL> avisc;
        
        //! ��ֵ�ճ�ȹ�������
        vector<OCP_DBL> bvisc;
};

/** @brief ����Lohrenz-Bray-Clark ճ�ȼ��㹫ʽ
*/
/// Lohrenz-Bray-Clark formula 
class ViscosityMethod03 : public ViscosityMethod
{
public:
    ViscosityMethod03(const ComponentParam& param, const USI& tarId);
    /// ����ճ��
    OCP_DBL CalViscosity(const ViscosityParams& vp) override;
    /// ����ճ�Ⱥ���Ӧ�ĵ���
    OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;

protected:
    USI             nc;       ///< �����
    vector<OCP_DBL> coef;     ///< LBCϵ��
    vector<OCP_DBL> Tc;       ///< ����ٽ��¶�
    vector<OCP_DBL> Pc;       ///< ����ٽ�ѹ��
    vector<OCP_DBL> Vcvis;    ///< ���ճ���ٽ����
    vector<OCP_DBL> MWC;      ///< ��ַ�������
    OCP_DBL MW;               ///< ���������
    vector<OCP_DBL> sqrtMWC;  ///< ��ַ���������ƽ����
    OCP_DBL xPc, xTc, xVc;    ///< ��������
    vector<OCP_DBL> auxA;     ///< ��������
    vector<OCP_DBL> auxB;     ///< ��������
};


/** @brief ճ�ȼ���ӿ���
*/
class ViscosityCalculation
{
public:
    /// Ĭ�Ϲ��캯��
    ViscosityCalculation() = default;
    /// ��װ
    void Setup(const ComponentParam& param, const USI& tarId);
    /// ����ճ��
    OCP_DBL CalViscosity(const ViscosityParams& vp) {
        return vM->CalViscosity(vp);
    }
    /// ����ճ�Ⱥ���Ӧ�ĵ���
    OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) {
        return vM->CalViscosity(vp, muP, muT, mux);
    }
protected:
    ViscosityMethod* vM;  ///< ճ�ȼ��㷽����
};


/** @brief �ʼ��㷽��
*/
class EnthalpyMethod
{
public:
    /// Ĭ�Ϲ��캯��
    EnthalpyMethod() = default;
    /// ������
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
    /// �����ʺ���Ӧ�ĵ���
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const = 0;
};


/** @brief Һ���ʺͼ��ʼ��㷽��
*/
class EnthalpyMethod01 : public EnthalpyMethod
{
public:
    /// ���캯��
    EnthalpyMethod01(const OCP_DBL& Trefin, const vector<OCP_DBL>& cpl1in, const vector<OCP_DBL>& cpl2in,
        const vector<OCP_DBL>& cpl3in, const vector<OCP_DBL>& cpl4in);
    /// ������
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
    /// �����ʺ���Ӧ�ĵ���
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
    OCP_DBL         Tref;  ///< �ο��¶�
    USI             nc;    ///< �����
    vector<OCP_DBL> cpl1;  ///< Һ���ʼ���ϵ��
    vector<OCP_DBL> cpl2;  ///< Һ���ʼ���ϵ��
    vector<OCP_DBL> cpl3;  ///< Һ���ʼ���ϵ��
    vector<OCP_DBL> cpl4;  ///< Һ���ʼ���ϵ��
};


/** @brief �����ʼ��㷽��
*/
class EnthalpyMethod02 : public EnthalpyMethod
{
public:
    /// ���캯��
    EnthalpyMethod02(const OCP_DBL& Trefin, const vector<OCP_DBL>& Tcritin,
        const vector<OCP_DBL>& cpg1in, const vector<OCP_DBL>& cpg2in,
        const vector<OCP_DBL>& cpg3in, const vector<OCP_DBL>& cpg4in,
        const vector<OCP_DBL>& hvaprin, const vector<OCP_DBL>& hvrin,
        const vector<OCP_DBL>& evin);
    /// ������
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
    /// �����ʺ���Ӧ�ĵ���
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
    OCP_DBL         Tref;    ///< �ο��¶�
    vector<OCP_DBL> Tcrit;   ///< ����ٽ��¶�
    USI             nc;      ///< �����
    vector<OCP_DBL> cpg1;    ///< �����ʼ���ϵ��
    vector<OCP_DBL> cpg2;    ///< �����ʼ���ϵ��
    vector<OCP_DBL> cpg3;    ///< �����ʼ���ϵ��
    vector<OCP_DBL> cpg4;    ///< �����ʼ���ϵ��
    vector<OCP_DBL> hvapr;   ///< �����ʼ���ϵ��
    vector<OCP_DBL> hvr;     ///< �����ʼ���ϵ��
    vector<OCP_DBL> ev;      ///< �����ʼ���ϵ��
};


/** @brief �ʼ���ӿ�
*/
class EnthalpyCalculation
{
public:
    /// ��װ����
    void Setup(const ComponentParam& param, const USI& tarId);
    /// ������
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const { return eM->CalEnthalpy(T, zi); }
    /// �����ʺ���Ӧ�ĵ���
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const { return eM->CalEnthalpy(T, zi, HT, Hz); }
protected:
    EnthalpyMethod* eM; ///< �ʼ��㷽����
};