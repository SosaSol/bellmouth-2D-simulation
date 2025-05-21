#!/usr/bin/env python3
"""
Post-process OpenFOAM outputs in parallel.:
- Extract head losses, mass flow rates, and max yPlus from each simulation case.
- Assemble results into CSV tables.
- Optionally plot iteration data.

Directory structure:
  ./outputs/<caseName>/postProcessing
Outputs written to:
  ./postProcessingOutputs/
"""

from postProcess import *
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union

start_time = time.time()

# ────────────────────────────────
# Main Workflow
# ────────────────────────────────

def main(write: bool = False, plot: bool = False) -> None:
    base_outputs = Path("outputs")
    if not base_outputs.exists():
        logger.error(f"Base outputs directory not found: {base_outputs}")
        sys.exit(1)

    all_subdirs = [d for d in base_outputs.iterdir() if d.is_dir()]
    if not all_subdirs:
        logger.error("No subdirectories found in outputs/.")
        sys.exit(1)

    summary_data = []

    for subdir in all_subdirs:
        output_dir = subdir
        postproc_dir = Path("postProcessingOutputs") / subdir.name

        cases = sorted([p for p in output_dir.iterdir() if p.is_dir()], key=natural_key)
        if not cases:
            logger.info("")
            logger.info(f"No cases found in {output_dir}")
            continue
        logger.info("")
        logger.info(f">>> Processing '{subdir.name}' with {len(cases)} cases")

        results: Dict[Tuple[str, str], Dict[Tuple[int, int], Dict[str, float]]] = {}

        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_single_case, case, postproc_dir, plot, write) for case in cases]
            for future in as_completed(futures):
                try:
                    result = future.result()
                except Exception:
                    logger.error("Worker crashed")
                    traceback.print_exc()
                    continue

                if not result:
                    continue

                (geom, key), dims, data = result

                logger.info("")
                logger.info(f"Processing: {data['case_name']}")
                logger.info(f"   - Iterations:      {data['Mb']}")
                logger.info(f"   - Pressure Out:    {data['head_loss']:.4f} Pa/m")
                logger.info(f"   - Head Loss:       {data['head_loss']:.4f} Pa/m")
                logger.info(f"   - Mass Flow Rate:  {data['massflow']:.4f} kg/s/m")
                logger.info(f"   - Max yPlus:       {data['max_yplus']:.2f}")
                if data['max_yplus'] > 5:
                    logger.warning("yPlus > 5!")
                elif data['max_yplus'] >= 1:
                    logger.warning("yPlus between 1 and 5.")

                results.setdefault((geom, key), {})[dims] = {
                    "head_loss": data["head_loss"],
                    "massflow":  data["massflow"],
                    "max_yplus": data["max_yplus"],
                }
                summary_data.append(data)

        if write and results:
            import pandas as pd
            postproc_dir.mkdir(parents=True, exist_ok=True)

            logger.info("")
            for (geom, key), data_dict in results.items():
                Mb_vals = sorted(set(mb for mb, _ in data_dict))
                Mw_vals = sorted(set(mw for _, mw in data_dict))
                idx = pd.Index(Mb_vals, name="Mb")
                cols = pd.Index(Mw_vals, name="Mw")

                df_head = pd.DataFrame(index=idx, columns=cols, dtype=float)
                df_flow = df_head.copy()
                df_yplus = df_head.copy()

                for (mb, mw), vals in data_dict.items():
                    df_head.at[mb, mw] = vals["head_loss"]
                    df_flow.at[mb, mw] = vals["massflow"]
                    df_yplus.at[mb, mw] = vals["max_yplus"]

                base = f"{geom}" if geom == "IN" else f"{geom}-{key}"
                df_head.to_csv(postproc_dir / f"{base}_head_losses.csv")
                df_flow.to_csv(postproc_dir / f"{base}_massflow.csv")
                df_yplus.to_csv(postproc_dir / f"{base}_max_yplus.csv")
                logger.info(f"Saved CSVs for {base} in {postproc_dir}")

        if write and summary_data:
            import pandas as pd
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(postproc_dir / "all_cases_summary.csv", index=False)
            logger.info("Saved all-cases summary to {postproc_dir}/all_cases_summary.csv")

    logger.info("")
    logger.info("All processing complete.")
    logger.info(f"Total runtime: {time.time() - start_time:.2f} seconds")

# ────────────────────────────────
# Script Entry Point
# ────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    main(write=args.write, plot=args.plot)