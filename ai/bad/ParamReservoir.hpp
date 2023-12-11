/*! \file    ParamReservoir.hpp
 *  \brief   �Ͳز����������ģ��
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMRESERVOIR_HEADER__
#define __PARAMRESERVOIR_HEADER__

 // Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/// һά����Ͳ������ݽṹ
class TableSet
{
public:
    /// ��ӡ�����Ļ
    void DisplayTable() const;

public:
    string                          name;       ///< �������
    USI                             colNum;     ///< ����
    vector<vector<vector<OCP_DBL>>> data;       ///< ����
};

/// ��ά����Ͳ������ݽṹ
class Table2
{
public:
    /// ���캯��
    Table2(const USI& n) { data.resize(n); }
    /// ���ñ�������
    void SetColNum() { colNum = data[0].size(); }
public:

    string                          refName;       ///< ��һά��������
    vector<OCP_DBL>                 refData;       ///< ��һά��������
    USI                             colNum;        ///< �ڶ�ά����
    vector<vector<vector<OCP_DBL>>> data;          ///< �������
};


/// ��ά�������
class Table2Set
{
public:
    string           name;   ///< �������
    vector<Table2>   data;   ///< �������
};

/// ����ʧ����
class HLoss
{
public:
    OCP_BOOL ifHLoss{ OCP_FALSE }; ///< �Ƿ�ʹ��
    OCP_BOOL obUse{ OCP_FALSE };   ///< �Ƿ�ʹ�ø�������ʧ
    OCP_DBL  obC{ -1 };            ///< ��������ʧ����
    OCP_DBL  obK{ -1 };            ///< ��������ʧ���� 
    OCP_BOOL ubUse{ OCP_FALSE };   ///< �Ƿ�ʹ�õͲ�����ʧ 
    OCP_DBL  ubC{ -1 };            ///< �ײ�����ʧ���� 
    OCP_DBL  ubK{ -1 };            ///< �ײ�����ʧ���� 
};

/// ���Բ���
class RockParam
{
public:
    string   type{ "LINEAR" };      ///< ѹ������
    OCP_DBL  Pref{ 14.7 };          ///< �ο�ѹ��
    OCP_DBL  Tref{ 60 };            ///< �ο��¶�
    OCP_DBL  cp1{ 3.406E-6 };       ///< һ��ѹ��ϵ��
    OCP_DBL  cp2{ 0 };              ///< ����ѹ��ϵ��
    OCP_DBL  ct{ 0 };               ///< ������ϵ��
    OCP_DBL  cpt{ 0 };              ///< ������ϵ��
    OCP_BOOL ConstRock{ OCP_TRUE }; ///< ��ʯ����Ƿ�Ϊ����
    OCP_DBL HCP1{ 35 };             ///< ��ʯ��ϵ��
    OCP_DBL HCP2{ 0 };              ///< ��ʯ��ϵ��
};


/// ��ʼ�Ͳ�ƽ�����
class EQUILParam
{
public:
    vector<OCP_DBL> data; ///< ƽ�����
};


/// Brooks-Corey �����͸��ë����ģ��
class BrooksCoreyParam
{
public:
    OCP_DBL sw_imm;     ///< �����ƶ�����ʪ�౥�Ͷ�
    OCP_DBL sn_imm;     ///< �����ƶ��ķ���ʪ�౥�Ͷ�
    OCP_DBL Pentry;     ///< �������ѹ��
    OCP_DBL Pcmax;      ///< �����ˮë����
    OCP_DBL Cw_kr;      ///< ��ʪ�������͸��ָ��
    OCP_DBL Cn_kr;      ///< ����ʪ�������͸��ָ��
    OCP_DBL C_pc;       ///< ë����ָ��
};


/// �߽���������
class BoundaryParam
{
public:
    /// ���캯��
    BoundaryParam(const string& Name) : name(Name) {}

    string   name;                     ///< �߽�����
    OCP_BOOL constP{ OCP_FALSE };      ///< �Ƿ�ѹ����
    OCP_DBL  P;                        ///< ��ѹ���Ƶ�ѹ��ֵ
};


/// �����ṹ1
template <typename T>
class Type_A_r
{
public:
    OCP_BOOL  activity{ OCP_FALSE }; ///< �Ƿ�ʹ��
    vector<T> data;                  ///< ��������
};

/// ��������
class ComponentParam
{
public:
    // Basic params
    /// ��ʼ������
    void Init();
    /// ���������Ϣ
    void InputCOMPONENTS(ifstream& ifs, const string& keyword);
    /// ��ֱ���ָ�������1
    Type_A_r<vector<OCP_DBL>>* FindPtr01(const string& varName);
    /// ����ο�ѹ���Ͳο��¶�
    void InputRefPR(ifstream& ifs, const string& keyword);
    /// ��ֱ���ָ�������2
    vector<OCP_DBL>* FindPtr02(const string& varName);
    /// �������������
    void InputCNAMES(ifstream& ifs);
    /// ����LBCϵ��
    void InputLBCCOEF(ifstream& ifs);
    /// ����BICϵ��
    void InputBIC(ifstream& ifs);
    /// �����ȶ��Է�����SSM����
    void InputSSMSTA(ifstream& ifs);
    /// �����ȶ��Է�����NR����
    void InputNRSTA(ifstream& ifs);
    /// ��������Ѽ����SSM����
    void InputSSMSP(ifstream& ifs);
    /// ��������Ѽ����NR����
    void InputNRSP(ifstream& ifs);
    /// ����RR���̵�������
    void InputRR(ifstream& ifs);

public:
    USI            NTPVT;               ///< PVT����
    USI            numCom{ 0 };         ///< ��ƽ�����������
    USI            numPhase{ 2 };       ///< ��ƽ�������������
    vector<string> Cname;               ///< ��ֵ�����
    Type_A_r<vector<OCP_DBL>> Tc;       ///< ����ٽ��¶�
    Type_A_r<vector<OCP_DBL>> Pc;       ///< ����ٽ�ѹ��
    Type_A_r<vector<OCP_DBL>> Vc;       ///< ����ٽ����
    Type_A_r<vector<OCP_DBL>> Zc;       ///< ����ٽ�ѹ������
    Type_A_r<vector<OCP_DBL>> MW;       ///< ��ַ�������
    Type_A_r<vector<OCP_DBL>> Acf;      ///< ���ƫ������
    Type_A_r<vector<OCP_DBL>> OmegaA;   ///< PR-EoSϵ��
    Type_A_r<vector<OCP_DBL>> OmegaB;   ///< PR-EoSϵ��
    Type_A_r<vector<OCP_DBL>> Vshift;   ///< �������任ϵ��
    Type_A_r<vector<OCP_DBL>> parachor; ///< parachor

    // for viscosity calculation
    Type_A_r<vector<OCP_DBL>> Vcvis;    ///< ճ�ȼ����ٽ����
    Type_A_r<vector<OCP_DBL>> Zcvis;    ///< ճ�ȼ���ѹ������
    vector<OCP_DBL>           LBCcoef;  ///< LBCճ�ȼ���ϵ��
    vector<vector<OCP_DBL>>   BIC;      ///< BIC����

    // Thermal only
    Type_A_r<vector<OCP_DBL>> molden;   ///< �ο��¶�ѹ���µ�Ħ��Ũ��
    Type_A_r<vector<OCP_DBL>> cp;       ///< ���ѹ��ϵ��
    Type_A_r<vector<OCP_DBL>> ct1;      ///< ��һ������ϵ��
    Type_A_r<vector<OCP_DBL>> ct2;      ///< ��һ������ϵ��
    Type_A_r<vector<OCP_DBL>> cpt;      ///< �¶�ѹ������ϵ��
    Type_A_r<vector<OCP_DBL>> cpl1;     ///< Һ���ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpl2;     ///< Һ���ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpl3;     ///< Һ���ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpl4;     ///< Һ���ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpg1;     ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpg2;     ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpg3;     ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> cpg4;     ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> hvapr;    ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> hvr;      ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> ev;       ///< �����ʼ���ϵ��
    Type_A_r<vector<OCP_DBL>> avisc;    ///< Һ����ֵ�ճ�ȹ������� 
    Type_A_r<vector<OCP_DBL>> bvisc;    ///< Һ����ֵ�ճ�ȹ������� 
    Type_A_r<vector<OCP_DBL>> avg;      ///< ������ֵ�ճ�ȹ������� 
    Type_A_r<vector<OCP_DBL>> bvg;      ///< ������ֵ�ճ�ȹ������� 


    Table2Set viscTab;                  ///< ճ���¶�ѹ��������

    vector<OCP_DBL> Pref;               ///< �ο�ѹ��
    vector<OCP_DBL> Tref;               ///< �ο��¶�
    vector<string> SSMparamSTA;         ///< �ȶ��Է���SSM����
    vector<string> NRparamSTA;          ///< �ȶ��Է���NR����
    vector<string> SSMparamSP;          ///< ����Ѽ���SSM����
    vector<string> NRparamSP;           ///< ����Ѽ���NR����
    vector<string> RRparam;             ///< RR������
};


/// ���ܿ��Ʋ���
class Miscstr
{
public:
    OCP_BOOL        ifMiscible{ OCP_FALSE };///< �Ƿ�ʹ�û���
    OCP_DBL         surTenRef{ -1 };        ///< �ο���������
    OCP_DBL         surTenEpt{ -1 };        ///< ����������ز���
    OCP_DBL         surTenPc{ -1 };         ///< ����������ز���
    OCP_DBL         surTenExp{ 0.25 };      ///< ����������ز���
};


/// �Ͳ��������ģ��
class ParamReservoir
{

public:

    string                   unitType;    ///< ��λ����
                                          
    OCP_DBL                  rsTemp;      ///< �Ͳس�ʼ�¶�
    vector<RockParam>        rockSet;     ///< ��ʯ������
    HLoss                    hLoss;       ///< ����ʧ������
    Miscstr                  miscstr;     ///< ���ܿ��Ʋ�����
    vector<BrooksCoreyParam> BCparam;     ///< Brooks-Coreyģ�Ͳ�����
    vector<BoundaryParam>    BDparam;     ///< �߽�����������

    Type_A_r<OCP_DBL> density;             ///< ��̬�¸����ܶ�
    Type_A_r<OCP_DBL> gravity;             ///< ��̬�¸�������ϵ��
    OCP_BOOL          ifThcon{ OCP_FALSE };///< �Ƿ�ʹ���ȴ���
    OCP_DBL           thcono{ 24 };        ///< �����ȴ���ϵ��
    OCP_DBL           thcong{ 24 };        ///< �����ȴ���ϵ��
    OCP_DBL           thconw{ 24 };        ///< ˮ���ȴ���ϵ��
    OCP_DBL           thconr{ 24 };        ///< ��ʯ�ȴ���ϵ��

    // mixture Models
    OCP_BOOL blackOil{ OCP_FALSE }; ///< �Ƿ�ʹ�ú���ģ��
    OCP_BOOL comps{ OCP_FALSE };    ///< �Ƿ�ʹ�����ģ��
    OCP_BOOL thermal{ OCP_FALSE };  ///< �Ƿ�ʹ������ģ��
    OCP_BOOL oil{ OCP_FALSE };      ///< �����Ƿ����
    OCP_BOOL gas{ OCP_FALSE };      ///< �����Ƿ����
    OCP_BOOL water{ OCP_FALSE };    ///< ˮ���Ƿ����
    OCP_BOOL disGas{ OCP_FALSE };   ///< �����ܽ����Ƿ����

    // flow model
    OCP_BOOL GRAVDR{ OCP_FALSE };   ///< �Ƿ�Ҫ��˫��ģ����ʹ����������

    ComponentParam comsParam;       ///< ��ֱ���������

    // SAT Region & PVT Region
    USI               NTSFUN{ 1 }; ///< ���Ͷ�������
    USI               NTPVT{ 1 };  ///< PVT����
    USI               NTROOC{ 1 }; ///< ��ʯ����

    TableSet SWFN_T;               ///< SWFN���
    TableSet SWOF_T;               ///< SWOF���
    TableSet SGFN_T;               ///< SGFN���
    TableSet SGOF_T;               ///< SGOF���
    TableSet SOF3_T;               ///< SOF3���
    TableSet PBVD_T;               ///< PBVD���
    // initial zi vs depth
    TableSet           ZMFVD_T;  ///< ZMFVD���
    TableSet           TEMPVD_T; ///< TEMPVD���
    vector<EQUILParam> EQUIL;    ///< ƽ����������

    // PVT properties
    USI numPhase;     ///< ����
    USI numCom;       ///< �����
    TableSet PVCO_T;  ///< PVCO���
    TableSet PVDO_T;  ///< PVDO���
    TableSet PVCDO_T; ///< PVCDO���
    TableSet PVDG_T;  ///< PVDG���
    TableSet PVTW_T;  ///< PVTW���


    OCP_BOOL  GARCIAW{ OCP_FALSE };   ///< �Ƿ�ʹ��GARCIAW����
    Table2Set PVTH2O;                 ///< PVTH2O���
    Table2Set PVTCO2;                 ///< PVTCO2���
    OCP_DBL   Psurf;                  ///< �ر�ѹ��
    OCP_DBL   Tsurf;                  ///< �ر��¶�

public:

    /// ���Ѱַ����1
    TableSet* FindPtrTable(const string& varName);
    /// ���Ѱַ����2
    Table2Set* FindPtrTable2(const string& varName);

    /// �����Ͳز�����ʼ��
    void Init();

    /// �Ͳر���ʼ��
    void InitTable();

    /// ������ֲ���
    void InputCOMPS(ifstream& ifs);

    /// �����Ͳس�ʼ�¶�
    void InputRTEMP(ifstream& ifs);

    /// ����һά���
    void InputTABLE(ifstream& ifs, const string& tabName);
    /// �����ά���
    void InputTABLE2(ifstream& ifs, const string& tabName);

    /// ���������ʯ����
    void InputROCK(ifstream& ifs);
    /// ����ǵ�����ʯ����
    void InputROCKT(ifstream& ifs);
    /// ��������ʧ����
    void InputHLOSS(ifstream& ifs);
    /// ����Brooks-Coreyģ�Ͳ���
    void InputBrooksCorey(ifstream& ifs);

    /// ������ܿ��Ʋ���
    void InputMISCSTR(ifstream& ifs);

    /// ������������
    void InputGRAVITY(ifstream& ifs);

    /// �����ܶȲ���
    void InputDENSITY(ifstream& ifs);

    /// �����ȴ�������
    void InputTHCON(ifstream& ifs, const string& keyword);

    /// ����ƽ���ʼ������
    void InputEQUIL(ifstream& ifs);

    /// ����������
    void InputTABDIMS(ifstream& ifs);

    /// ���������
    void InputNCOMPS(ifstream& ifs) {
        vector<string> vbuf;
        ReadLine(ifs, vbuf);
        comsParam.numCom = stoi(vbuf[0]);
        numCom = comsParam.numCom;
    }
    /// �����������
    void InputCNAMES(ifstream& ifs) { comsParam.InputCNAMES(ifs); };
    /// ������ֲ���
    void InputCOMPONENTS(ifstream& ifs, const string& keyword)
    {
        comsParam.InputCOMPONENTS(ifs, keyword);
    }
    /// ����LBCճ�Ȳ���
    void InputLBCCOEF(ifstream& ifs) { comsParam.InputLBCCOEF(ifs); }
    /// ����BIC����
    void InputBIC(ifstream& ifs) { comsParam.InputBIC(ifs); };
    /// ������ֲο�ѹ��
    void InputRefPR(ifstream& ifs, const string& keyword)
    {
        comsParam.InputRefPR(ifs, keyword);
    };

    /// �����ȶ��Է���SSM����
    void InputSSMSTA(ifstream& ifs) { comsParam.InputSSMSTA(ifs); };
    /// �����ȶ��Է���NR����
    void InputNRSTA(ifstream& ifs) { comsParam.InputNRSTA(ifs); };
    /// ��������Ѽ���SSM����
    void InputSSMSP(ifstream& ifs) { comsParam.InputSSMSP(ifs); };
    /// ��������Ѽ���NR����
    void InputNRSP(ifstream& ifs) { comsParam.InputNRSP(ifs); };
    /// ����RR����������
    void InputRR(ifstream& ifs) { comsParam.InputRR(ifs); };

    /// ����߽����
    void InputBoundary(ifstream& ifs);

    /// ������������ȷ��
    void CheckParam();

    /// �����ʯ����������ȷ��
    void CheckRock();

    /// ���CPL����������ȷ��
    void CheckCPL();

    /// ���CPL����������ȷ��
    void CheckCPG();
};

#endif /* end if __PARAMRESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/