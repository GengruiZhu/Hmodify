
import os
import subprocess
import configparser
import logging
from pathlib import Path

# 设置日志，记录到总输出目录的 modify_hic.log
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_command(cmd, error_msg):
    """运行 shell 命令，失败时记录错误并退出"""
    try:
        result = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
        logger.info(f"命令执行成功: {cmd}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"{error_msg}: {e.stderr}")
        exit(1)

def extract_chromosome(agp_file, start_utg, end_utg, output_file):
    """提取 AGP 文件中从 start_utg 到 end_utg 的染色体记录"""
    logger.info(f"提取染色体从 {start_utg} 到 {end_utg} 到 {output_file}")
    cmd = f"""awk -v start="{start_utg}" -v end="{end_utg}" '
        $1 == start || $6 == start {{flag=1}}
        flag {{print; if ($1 == end || $6 == end) exit}}' "{agp_file}" > "{output_file}" """
    run_command(cmd, f"提取染色体从 {start_utg} 到 {end_utg} 失败")
    if not os.path.getsize(output_file):
        logger.error(f"错误：{output_file} 为空")
        exit(1)

def extract_utg_range(agp_file, start_utg, end_utg, output_file):
    """提取 AGP 文件中从 start_utg 到 end_utg 的复制片段"""
    logger.info(f"提取复制片段 {start_utg}-{end_utg} 到 {output_file}")
    cmd = f"""awk -v start="{start_utg}" -v end="{end_utg}" '
        $1 == start || $6 == start {{flag=1}}
        flag {{print; if ($1 == end || $6 == end) exit}}' "{agp_file}" > "{output_file}" """
    run_command(cmd, f"提取片段 {start_utg}-{end_utg} 失败")
    if not os.path.getsize(output_file):
        logger.error(f"错误：{output_file} 为空")
        exit(1)

def extract_reference_utg(agp_file, ref_utg, output_file):
    """提取参考染色体的 utg 记录"""
    logger.info(f"提取参考 utg {ref_utg} 到 {output_file}")
    cmd = f"awk -v utg='{ref_utg}' '$1 == utg || $6 == utg {{print}}' '{agp_file}' > '{output_file}'"
    run_command(cmd, f"提取参考 utg {ref_utg} 失败")
    if not os.path.getsize(output_file):
        logger.error(f"错误：{output_file} 为空")
        exit(1)

def insert_utg_segment(chr_file, insert_utg, fragment_file, output_file):
    """在指定 utg 后插入复制片段，保留后续记录"""
    logger.info(f"在 {chr_file} 的 {insert_utg} 后插入 {fragment_file} 到 {output_file}")
    cmd = f'''
    awk -F'\\t' -v insert_utg="{insert_utg}" '
        BEGIN {{output=1}}
        {{if (output) print > "{output_file}.tmp1"; else print > "{output_file}.tmp2"}}
        $1 == insert_utg || $6 == insert_utg {{output=0}}
    ' "{chr_file}"
    '''
    run_command(cmd, f"分割 {chr_file} 失败")
    # 合并：前半部分 + 插入片段 + 后半部分
    cmd = f"""
    cat "{output_file}.tmp1" "{fragment_file}" "{output_file}.tmp2" > "{output_file}"
    rm -f "{output_file}.tmp1" "{output_file}.tmp2"
    """
    run_command(cmd, f"合并到 {output_file} 失败")
    if not os.path.getsize(output_file):
        logger.error(f"错误：{output_file} 为空")
        exit(1)

def main():
    # 读取配置文件
    logger.info("读取 config.txt")
    config = configparser.ConfigParser()
    config.read("config.txt")

    # 获取全局参数
    output_dir = config["DEFAULT"].get("OUTPUT_DIR", "")
    agp_file = config["DEFAULT"].get("AGP_FILE", "")
    fasta_file = config["DEFAULT"].get("FASTA_FILE", "")

    # 验证全局参数
    if not all([output_dir, agp_file, fasta_file]):
        logger.error("错误：缺少全局参数（OUTPUT_DIR, AGP_FILE, FASTA_FILE）")
        exit(1)

    # 创建总输出目录并设置日志
    os.makedirs(output_dir, exist_ok=True)
    handler = logging.FileHandler(f"{output_dir}/modify_hic.log", mode="a")
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] - %(message)s"))
    logger.addHandler(handler)

    # 检查输入文件
    logger.info("检查输入文件")
    for file in [agp_file, fasta_file]:
        if not os.path.exists(file) or not os.path.isfile(file) or not os.access(file, os.R_OK):
            logger.error(f"错误：无法访问文件 {file}")
            exit(1)
        logger.info(f"{file}: {os.path.getsize(file)} 字节")

    # 处理每组操作
    for section in config.sections():
        logger.info(f"处理操作组：{section}")
        params = {
            "CHR_A": config[section].get("CHR_A", ""),
            "START_UTG_FULL_A": config[section].get("START_UTG_FULL_A", ""),
            "END_UTG_FULL_A": config[section].get("END_UTG_FULL_A", ""),
            "START_UTG_A": config[section].get("START_UTG_A", ""),
            "END_UTG_A": config[section].get("END_UTG_A", ""),
            "INSERT_AFTER_UTG_B": config[section].get("INSERT_AFTER_UTG_B", ""),
            "CHR_B": config[section].get("CHR_B", ""),
            "START_UTG_FULL_B": config[section].get("START_UTG_FULL_B", ""),
            "END_UTG_FULL_B": config[section].get("END_UTG_FULL_B", ""),
            "START_UTG_B": config[section].get("START_UTG_B", ""),
            "END_UTG_B": config[section].get("END_UTG_B", ""),
            "INSERT_AFTER_UTG_A": config[section].get("INSERT_AFTER_UTG_A", ""),
            "REF_CHR": config[section].get("REF_CHR", ""),
            "REF_UTG": config[section].get("REF_UTG", "")
        }

        # 验证必填参数
        if not all([params["CHR_A"], params["CHR_B"], params["START_UTG_FULL_A"], params["END_UTG_FULL_A"],
                    params["START_UTG_FULL_B"], params["END_UTG_FULL_B"]]):
            logger.error(f"错误：{section} 缺少 CHR_A, CHR_B, START_UTG_FULL_A, END_UTG_FULL_A, START_UTG_FULL_B, END_UTG_FULL_B")
            continue
        if not ((params["START_UTG_A"] and params["END_UTG_A"] and params["INSERT_AFTER_UTG_B"]) or
                (params["START_UTG_B"] and params["END_UTG_B"] and params["INSERT_AFTER_UTG_A"])):
            logger.error(f"错误：{section} 至少需要一组复制参数（A->B 或 B->A）")
            continue

        # 创建 PartXX 子目录
        part_dir = f"{output_dir}/{section}"
        os.makedirs(part_dir, exist_ok=True)
        logger.info(f"Part 目录: {part_dir}")

        # 检查 AGP 中的 utg
        utgs = [u for u in [params["START_UTG_FULL_A"], params["END_UTG_FULL_A"], params["START_UTG_A"],
                            params["END_UTG_A"], params["INSERT_AFTER_UTG_B"], params["START_UTG_FULL_B"],
                            params["END_UTG_FULL_B"], params["START_UTG_B"], params["END_UTG_B"],
                            params["INSERT_AFTER_UTG_A"], params["REF_UTG"]] if u]
        for utg in utgs:
            cmd = f"grep -q '{utg}' '{agp_file}'"
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError:
                logger.error(f"错误：在 AGP 文件中未找到 UTG {utg}")
                exit(1)

        # 步骤1: 提取完整染色体
        extract_chromosome(agp_file, params["START_UTG_FULL_A"], params["END_UTG_FULL_A"],
                           f"{part_dir}/Chromosome{params['CHR_A']}_full.txt")
        extract_chromosome(agp_file, params["START_UTG_FULL_B"], params["END_UTG_FULL_B"],
                           f"{part_dir}/Chromosome{params['CHR_B']}_full.txt")

        # 步骤2: 提取复制片段
        if params["START_UTG_A"] and params["END_UTG_A"]:
            extract_utg_range(agp_file, params["START_UTG_A"], params["END_UTG_A"],
                              f"{part_dir}/Chromosome{params['CHR_A']}_add.txt")
        if params["START_UTG_B"] and params["END_UTG_B"]:
            extract_utg_range(agp_file, params["START_UTG_B"], params["END_UTG_B"],
                              f"{part_dir}/Chromosome{params['CHR_B']}_add.txt")

        # 步骤3: 插入片段生成新染色体
        if params["START_UTG_A"] and params["END_UTG_A"] and params["INSERT_AFTER_UTG_B"]:
            logger.info(f"在 Chromosome{params['CHR_B']}_full 的 {params['INSERT_AFTER_UTG_B']} 后插入 Chromosome{params['CHR_A']}_add")
            insert_utg_segment(
                f"{part_dir}/Chromosome{params['CHR_B']}_full.txt",
                params["INSERT_AFTER_UTG_B"],
                f"{part_dir}/Chromosome{params['CHR_A']}_add.txt",
                f"{part_dir}/Chromosome{params['CHR_B']}_new.txt"
            )
        else:
            logger.info(f"不修改 Chromosome{params['CHR_B']}，保持原样")
            cmd = f"cp '{part_dir}/Chromosome{params['CHR_B']}_full.txt' '{part_dir}/Chromosome{params['CHR_B']}_new.txt'"
            run_command(cmd, f"复制 Chromosome{params['CHR_B']} 失败")

        if params["START_UTG_B"] and params["END_UTG_B"] and params["INSERT_AFTER_UTG_A"]:
            logger.info(f"在 Chromosome{params['CHR_A']}_full 的 {params['INSERT_AFTER_UTG_A']} 后插入 Chromosome{params['CHR_B']}_add")
            insert_utg_segment(
                f"{part_dir}/Chromosome{params['CHR_A']}_full.txt",
                params["INSERT_AFTER_UTG_A"],
                f"{part_dir}/Chromosome{params['CHR_B']}_add.txt",
                f"{part_dir}/Chromosome{params['CHR_A']}_new.txt"
            )
        else:
            logger.info(f"不修改 Chromosome{params['CHR_A']}，保持原样")
            cmd = f"cp '{part_dir}/Chromosome{params['CHR_A']}_full.txt' '{part_dir}/Chromosome{params['CHR_A']}_new.txt'"
            run_command(cmd, f"复制 Chromosome{params['CHR_A']} 失败")

        # 步骤4: 提取参考染色体（如果提供）
        if params["REF_UTG"] and params["REF_CHR"]:
            extract_reference_utg(agp_file, params["REF_UTG"],
                                 f"{part_dir}/Chromosome{params['REF_CHR']}.txt")

        # 步骤5: 组合染色体生成最终utg列表
        combined_file = f"{part_dir}/Chromosome"
        fasta_name = "Chromosome"
        components = []

        if params["REF_CHR"]:
            components.append(f"{part_dir}/Chromosome{params['REF_CHR']}.txt")
            combined_file += params["REF_CHR"]
            fasta_name += params["REF_CHR"]
        if params["START_UTG_A"] and params["END_UTG_A"] and params["INSERT_AFTER_UTG_B"]:
            components.append(f"{part_dir}/Chromosome{params['CHR_B']}_new.txt")
            combined_file += f"-{params['CHR_B']}new"
            fasta_name += f"-{params['CHR_B']}"
        else:
            components.append(f"{part_dir}/Chromosome{params['CHR_B']}_new.txt")
            combined_file += f"-{params['CHR_B']}"
            fasta_name += f"-{params['CHR_B']}"
        if params["START_UTG_B"] and params["END_UTG_B"] and params["INSERT_AFTER_UTG_A"]:
            components.append(f"{part_dir}/Chromosome{params['CHR_A']}_new.txt")
            combined_file += f"-{params['CHR_A']}new"
            fasta_name += f"-{params['CHR_A']}"
        elif not (params["START_UTG_A"] and params["END_UTG_A"] and params["INSERT_AFTER_UTG_B"]):
            components.append(f"{part_dir}/Chromosome{params['CHR_A']}_new.txt")
            combined_file += f"-{params['CHR_A']}"
            fasta_name += f"-{params['CHR_A']}"

        combined_file += ".txt"
        fasta_name += "_new.fasta"

        logger.info(f"组合染色体到 {combined_file}")
        cmd = f"cat {' '.join(components)} > '{combined_file}'"
        run_command(cmd, "生成最终utg列表失败")

        # 步骤6: 提取fasta序列
        logger.info(f"提取fasta序列到 {fasta_name}")
        # 修改：调试第6列内容，宽松正则，排除空行
        cmd = f"cut -f6 '{combined_file}' > '{part_dir}/raw_utg_patterns.txt'"
        run_command(cmd, "生成原始utg列表失败")
        cmd = f"cut -f6 '{combined_file}' | grep -E '^utg[0-9]+l$' | sort | uniq > '{part_dir}/utg_patterns.txt'"
        run_command(cmd, "生成utg_patterns.txt 失败")
        pattern_count = int(run_command(f"wc -l < '{part_dir}/utg_patterns.txt'", "统计模式数失败").strip())
        logger.info(f"找到 {pattern_count} 个唯一UTG模式")
        if pattern_count == 0:
            logger.error(f"错误：utg_patterns.txt 为空，检查 {part_dir}/raw_utg_patterns.txt")
            with open(f"{part_dir}/raw_utg_patterns.txt") as f:
                logger.error(f"原始第6列内容：\n{f.read()}")
            exit(1)
        cmd = f"seqkit grep -f '{part_dir}/utg_patterns.txt' '{fasta_file}' -o '{part_dir}/{fasta_name}'"
        run_command(cmd, f"生成 FASTA 失败")

        # 检查 FASTA 文件
        logger.info(f"检查 FASTA 文件: {fasta_name}")
        cmd = f"ls -lh '{part_dir}/{fasta_name}'"
        run_command(cmd, "显示 FASTA 文件失败")
        cmd = f"grep -c '^>' '{part_dir}/{fasta_name}'"
        sequence_count = int(run_command(cmd, "统计序列数失败").strip())
        logger.info(f"FASTA 文件中找到 {sequence_count} 条序列")
        if sequence_count != pattern_count:
            logger.warning(f"警告：UTG 模式数 ({pattern_count}) 与序列数 ({sequence_count}) 不一致")
            cmd = f"seqkit seq --name '{fasta_file}' | grep -v -f '{part_dir}/utg_patterns.txt' > '{part_dir}/missing_utgs.txt'"
            run_command(cmd, f"生成缺失utg记录失败")
            with open(f"{part_dir}/missing_utgs.txt") as f:
                logger.warning(f"缺失的 UTG：\n{f.read()}")

        logger.info(f"处理完成：{part_dir}/{fasta_name}")

    logger.removeHandler(handler)

if __name__ == "__main__":
    main()
