#!/usr/bin/env python3
"""
RIBLL2026 模拟主入口
用法:
  python3 simulation.py            # 启动GUI
  python3 simulation.py -n 10000   # 直接模拟10000个事件
"""

import sys
import os
from datetime import datetime

# 确保项目根目录在Python路径中
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import tomllib
except ImportError:
    import tomli as tomllib

from ribll_sim.simulation.engine import run_simulation, save_to_root


def load_config(config_path='config.toml'):
    with open(config_path, 'rb') as f:
        return tomllib.load(f)


def main():
    use_gui = True
    n_events_override = None

    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == '-n' and i + 1 < len(args):
            use_gui = False
            n_events_override = int(args[i + 1])
            i += 2
        elif args[i] == '-c' and i + 1 < len(args):
            config_path = args[i + 1]
            i += 2
        else:
            i += 1

    config = load_config()

    if n_events_override is not None:
        config.setdefault('simulation', {})['n_events'] = n_events_override

    if use_gui:
        print("启动GUI模式...")
        from ribll_sim.gui.app import launch_gui
        launch_gui(config)
    else:
        print("启动批量模拟模式...")
        results = run_simulation(config)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = os.path.join(
            config.get('output', {}).get('data_dir', './output'),
            timestamp
        )
        save_to_root(results, output_dir, timestamp)

        with open(os.path.join(output_dir, 'config.toml'), 'w') as f:
            import toml
            toml.dump(config, f)


if __name__ == '__main__':
    main()