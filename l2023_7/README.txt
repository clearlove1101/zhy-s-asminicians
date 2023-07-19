{
    Equations{
        EOS{
            Tillotson等状态方程
        },
        solid{
            ROCK等描述固体的方程
        },
        others{
            其他方程
        }
    },
    integrator{
        step{
            积分器，例如 v[t+1]= v[t] + a[t] * dt
        }
    },
    tools{
        tools{
            创建粒子时若需额外写函数，放在这里
        }
    },
    visualization{
        visual{
            可视化工具
        }
    },
    main{
        最终运行的文件
    }
}











