/*! \file    WellPeaceman.hpp 
 *  \brief   WellPeacema class declaration 
 *  \author  Shizhe Li 
 *  \date    Aug/17/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLPEACEMAN_HEADER__
#define __WELLPEACEMAN_HEADER__

// Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "Well.hpp"
#include "WellPerf.hpp"

using namespace std;

/// Peaceman��ģ��ģ��
PeacemanWell : public Well{
public:
    /// ���������Ϣ
    void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) override;

    /// ��װ��
    void Setup(const Bulk& bk, const vector<SolventINJ>& sols) override;
    
    /// ��ʼ����ѹ��
    void InitWellP(const Bulk& bk) override;

    /// ��龮����ģʽ
    void CheckOptMode(const Bulk& bk) override;

    /// ʱ�䲽��ʼ���㾮����
    void CalFluxInit(const Bulk& bk) override;

    /// ���㾮����
    void CalFlux(const Bulk& bk) override;

    /// ��龮ѹ��
    ReservoirState CheckP(const Bulk& bk) override;

    /// �������ע����
    void CalIPRate(const Bulk& bk, const OCP_DBL& dt) override;

    /// ����ʱ�䲽�����仯
    OCP_DBL CalMaxChangeTime() const override;

    /// ����ţ�ٲ������仯
    OCP_DBL CalMaxChangeNR() override;

    /// ����Ϊ��һ��ʱ�䲽��״̬
    void ResetToLastTimeStep(const Bulk& bk) override;

    /// ���浱ǰʱ�䲽��״̬
    void UpdateLastTimeStep() override;

protected:
    /// ���㾮ָ��
    void CalWI(const Bulk& bk);

    /// ���㾮������
    void CalTrans(const Bulk& bk);

    /// ���㾮����
    void CalFlux(const Bulk& bk, const OCP_BOOL ReCalXi);

    /// �������ע������
    OCP_DBL CalInjRateMaxBHP(const Bulk& bk);

    /// �������ע������
    OCP_DBL CalProdRateMinBHP(const Bulk& bk);

    /// ����ע������
    void CalInjQj(const Bulk& bk, const OCP_DBL& dt);

    /// ������������
    void CalProdQj(const Bulk& bk, const OCP_DBL& dt);

    /// ��齻����
    ReservoirState CheckCrossFlow(const Bulk& bk);

    /// ��������ע������
    void CalFactor(const Bulk& bk) const;

    /// �������ѹ��
    void CaldG(const Bulk& bk);

    /// ����ע�뾮���ѹ��
    void CalInjdG(const Bulk& bk);

    /// �������������ѹ��
    void CalProddG(const Bulk& bk);

    /// �������������ѹ���1
    void CalProddG01(const Bulk& bk);

    /// �������������ѹ���2
    void CalProddG02(const Bulk& bk);

    /// �������ѹ��
    void CalPerfP() { for (USI p = 0; p < numPerf; p++) perf[p].P = bhp + dG[p]; }

protected:

    vector<OCP_DBL> dG;               ///< ���ѹ��
    mutable vector<OCP_DBL> factor;   ///< ע����������
};

class PeacemanWellIsoT : public PeacemanWell
{
public:
    /// FIM�������㾮�в�
    void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const override;
    /// FIM������ȡ����
    void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// FIM����װ�侮����
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// FIM����װ��ע�뾮����
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// FIM����װ������������
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    /// IMPEC������ȡ����
    void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// IMPEC����װ�侮����
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// IMPEC����װ��ע�뾮����
    void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// IMPEC����װ������������
    void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

};


class PeacemanWellT : public PeacemanWell
{
public:
    /// FIM���������������в�
    void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const override;
    /// FIM������ȡ��������
    void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// FIM����װ������������
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// FIM����װ������ע�뾮����
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// FIM����װ����������������
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    /// ������
    void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) override { OCP_ABORT("NOT USED!"); }
    /// ������
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override { OCP_ABORT("NOT USED!"); }
};

#endif /* end if __WELLPEACEMAN_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/