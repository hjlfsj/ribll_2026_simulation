# 标准的 A(B,D)C型反应 的运动学曲线绘制
# 可以选择随机生成，或者指定重离子的质心系角度生成

#对照Lise++已经进行的测试：
# 14O(2H,2H)14O, 15MeV
# 14O(2H,2H)14O, 0MeV
# 14O(1H,1H)14O, 15MeV
# 14O(1H,1H)14O, 0MeV
# 14O(2H, 6Li)10C, 0MeV












import ROOT
import numpy as np
import time
from datetime import datetime
from array import array
from math import pi, sqrt, cos, sin, acos
import sys
sys.path.append('/home/ribll2026/ribll2026_www/github_code/ribll_2026_simulation/reference/pyroot/standard_tools')  # Linux/Mac



from generate_kinematics import get_reaction_parameters
from generate_kinematics import simulate_single_event



E_k = 35*14.00421
excit_e = 0

particle_A = '14O'
particle_B = '2H'
particle_C = '10C'
particle_D = '6Li'
# A (B,D)C，A, C设置为重离子

def main():
    ROOT.gROOT.SetBatch(True)
    base_params = get_reaction_parameters(
        projectile_mass = particle_A,
        projectile_energy = E_k,
        target_type = particle_B,
        ejectile1_mass = particle_C,
        ejectile1_excitation = excit_e,
        ejectile2_mass = particle_D
    )
    print(f"反应参数:")
    print(f"  入射粒子: {base_params['proj_label']}, E = {base_params['E_proj_kin']} MeV")
    print(f"  靶核: {base_params['target_label']}")
    print(f"  出射粒子1: {base_params['ej1_label']}, Ex = {base_params['excitation1']} MeV")
    print(f"  出射粒子2: {base_params['ej2_label']}")
    print(f"  Q值: {base_params['Q_value']:.3f} MeV")
# 
    start_time = time.time()
    
    # 存储所有角度的结果
    all_results = []
    
    theta_cm_min = 0.1
    theta_cm_max = 180
    n_valid = int((theta_cm_max - theta_cm_min)/0.1)
    # 循环所有质心系角度
    for i in range(n_valid):
        
        if i % 1000 == 0:
            elapsed = time.time() - start_time
            print(f"已完成 {i} 个事件，耗时 {elapsed:.2f} 秒")
           
        base_params = get_reaction_parameters(
        projectile_mass = particle_A,
        projectile_energy = E_k,
        target_type = particle_B,
        ejectile1_mass = particle_C,
        ejectile1_excitation = excit_e,
        ejectile2_mass = particle_D,
        theta_cm = theta_cm_min + i*0.1
            )
        # 计算运动学
        result = simulate_single_event(base_params)
        
        if result is not None:
            all_results.append(result)
    
    total_time = time.time() - start_time
    n_valid = len(all_results)
    
    print(f"\n计算完成!")
    print(f"总时间: {total_time:.2f} 秒")
    # 准备数据用于TGraph
    n_points = n_valid
    
    # 提取数据
    theta_ej1_lab_list = [r[0] for r in all_results]
    E_ej1_lab_list = [r[1] for r in all_results]
    theta_ej2_lab_list = [r[2] for r in all_results]
    E_ej2_lab_list = [r[3] for r in all_results]
    
    # 创建TGraph对象
    # 1. 出射粒子1实验室系能量-角度
    gr1 = ROOT.TGraph(n_points)
    gr1.SetTitle(f"{base_params['ej1_label']}* : E vs theta_lab  ;theta_lab [deg];E [MeV]")
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(0.5)
    gr1.SetMarkerColor(ROOT.kRed)
    gr1.SetLineColor(ROOT.kRed)
    gr1.SetLineWidth(2)
    
    # 2. 出射粒子2实验室系能量-角度
    gr2 = ROOT.TGraph(n_points)
    gr2.SetTitle(f"{base_params['ej2_label']} : E vs theta_lab  ;theta_lab [deg];E [MeV]")
    gr2.SetMarkerStyle(21)
    gr2.SetMarkerSize(0.5)
    gr2.SetMarkerColor(ROOT.kBlue)
    gr2.SetLineColor(ROOT.kBlue)
    gr2.SetLineWidth(2)
    
    # 3. 双粒子角度关联
    gr3 = ROOT.TGraph(n_points)
    gr3.SetTitle(f"Angular Correlation;theta_lab({base_params['ej1_label']}*) [deg];theta_lab({base_params['ej2_label']}) [deg]")
    gr3.SetMarkerStyle(22)
    gr3.SetMarkerSize(0.5)
    gr3.SetMarkerColor(ROOT.kGreen+2)
    gr3.SetLineColor(ROOT.kGreen+2)
    gr3.SetLineWidth(2)
    
    # 4. 出射粒子1质心系角度-实验室系能量
    gr4 = ROOT.TGraph(n_points)
    gr4.SetTitle(f"{base_params['ej1_label']}* CM Angle vs Lab Energy;theta_CM [deg];E_lab [MeV]")
    gr4.SetMarkerStyle(23)
    gr4.SetMarkerSize(0.5)
    gr4.SetMarkerColor(ROOT.kMagenta)
    gr4.SetLineColor(ROOT.kMagenta)
    gr4.SetLineWidth(2)
    
    # 5. 双粒子能量关联
    gr5 = ROOT.TGraph(n_points)
    gr5.SetTitle(f"Energy Correlation;E_lab({base_params['ej1_label']}*) [MeV];E_lab({base_params['ej2_label']}) [MeV]")
    gr5.SetMarkerStyle(24)
    gr5.SetMarkerSize(0.5)
    gr5.SetMarkerColor(ROOT.kOrange+1)
    gr5.SetLineColor(ROOT.kOrange+1)
    gr5.SetLineWidth(2)
    
    # 6. 实验室角度 vs 质心角度
    gr6 = ROOT.TGraph(n_points)
    gr6.SetTitle(f"Lab Angle vs CM Angle;theta_CM [deg];theta_lab [deg]")
    gr6.SetMarkerStyle(25)
    gr6.SetMarkerSize(0.5)
    gr6.SetMarkerColor(ROOT.kCyan+1)
    gr6.SetLineColor(ROOT.kCyan+1)
    gr6.SetLineWidth(2)
    
    # 填充TGraph数据
    for i, result in enumerate(all_results):
        (theta_ej1_lab, E_ej1_lab,
         theta_ej2_lab, E_ej2_lab,
         theta_ej1_cm, E_ej1_cm) = result
        
        gr1.SetPoint(i, theta_ej1_lab, E_ej1_lab)
        gr2.SetPoint(i, theta_ej2_lab, E_ej2_lab)
        gr3.SetPoint(i, theta_ej1_lab, theta_ej2_lab)
        gr4.SetPoint(i, theta_ej1_cm, E_ej1_lab)
        gr5.SetPoint(i, E_ej2_lab, E_ej1_lab)
        gr6.SetPoint(i, theta_ej1_cm, theta_ej1_lab)
    
    # 创建ROOT画布并绘制
    canvas = ROOT.TCanvas("canvas", "Kinematic Curves", 1600, 1200)
    canvas.Divide(3, 2)
    
    # 设置全局样式
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)
    
    # 绘制图1
    canvas.cd(1)
    gr1.Draw("AP")
    canvas.cd(1).Update()
    
    # 绘制图2
    canvas.cd(2)
    gr2.Draw("AP")
    canvas.cd(2).Update()
    
    # 绘制图3
    canvas.cd(3)
    gr3.Draw("AP")
    canvas.cd(3).Update()
    
    # 绘制图4
    canvas.cd(4)
    gr4.Draw("AP")
    canvas.cd(4).Update()
    
    # 绘制图5
    canvas.cd(5)
    gr5.Draw("AP")
    canvas.cd(5).Update()
    
    # 绘制图6
    canvas.cd(6)
    gr6.Draw("AP")
    canvas.cd(6).Update()
    
    canvas.Update()
    path_pdf = f"../data_pdf/"
    name_pdf = f"kinematic_curves_{particle_A}_{particle_B}_{particle_D}_{particle_C}_{excit_e:.1f}MeV.pdf"
    canvas.SaveAs(path_pdf + name_pdf)
    
    path_root = f"../data_root/"
    name_root = f"kinematic_curves_{particle_A}_{particle_B}_{particle_D}_{particle_C}_{excit_e:.1f}MeV.root"
    
    # 新建root 文件
    f = ROOT.TFile(path_root + name_root, "RECREATE")

    # 将TGraph写入文件
    gr1.Write("gr1")
    gr2.Write("gr2")
    gr3.Write("gr3")
    gr4.Write("gr4")
    gr5.Write("gr5")
    gr6.Write("gr6")

    # 关闭文件
    f.Close()
    print(f"\n已保存为kinematic_curves.root")
    
    
    
    
 
    print("\n" + "="*70)
    print("运动学扫描统计:")
    print("="*70)
    
    if n_valid > 0:
        # 计算统计量
        theta_ej1_lab_arr = np.array(theta_ej1_lab_list)
        E_ej1_lab_arr = np.array(E_ej1_lab_list)
        theta_ej2_lab_arr = np.array(theta_ej2_lab_list)
        E_ej2_lab_arr = np.array(E_ej2_lab_list)
        
        print(f"\n{base_params['ej1_label']}* 实验室系:")
        print(f"  角度范围: {theta_ej1_lab_arr.min():.2f} 到 {theta_ej1_lab_arr.max():.2f} 度")
        print(f"  能量范围: {E_ej1_lab_arr.min():.2f} 到 {E_ej1_lab_arr.max():.2f} MeV")
        
        print(f"\n{base_params['ej2_label']} 实验室系:")
        print(f"  角度范围: {theta_ej2_lab_arr.min():.1f} 到 {theta_ej2_lab_arr.max():.1f} 度")
        print(f"  能量范围: {E_ej2_lab_arr.min():.1f} 到 {E_ej2_lab_arr.max():.1f} MeV")
        
    
    print("\n" + "="*70)
    print("程序完成")

# ====================== 直接运行 ======================
if __name__ == "__main__":
    main()