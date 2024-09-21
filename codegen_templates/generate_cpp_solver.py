from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import json
import argparse
import os
import inspect

def gen_cmake(env, solver_data, solver_path):
    template = env.get_template("CMakeLists.txt")
    content = template.render(solver_data)

    path = Path(solver_path, "CMakeLists.txt")
    with open(str(path), mode="w", encoding="utf-8") as cmakelists:
        cmakelists.write(content)
        #print(f"{content}")

def gen_HS_cpp(env, solver_data, solver_path):
    template = env.get_template("HomotopySolver.cpp")
    content = template.render(solver_data)

    path = Path(solver_path, "HomotopySolver.cpp")
    with open(str(path), mode="w", encoding="utf-8") as HS_cpp:
        HS_cpp.write(content)
        #print(f"{content}")

def gen_HS_hpp(env, solver_data, solver_path):
    template = env.get_template("HomotopySolver.hpp")
    content = template.render(solver_data)

    path = Path(solver_path, "HomotopySolver.hpp")
    with open(str(path), mode="w", encoding="utf-8") as HS_hpp:
        HS_hpp.write(content)
        #print(f"{content}")

def gen_PO_hpp(env, problem_options, solver_path):
    template = env.get_template("ProblemOptions.hpp")
    content = template.render(problem_options)

    path = Path(solver_path, "ProblemOptions.hpp")
    with open(str(path), mode="w", encoding="utf-8") as PO_hpp:
        PO_hpp.write(content)
        #print(f"{content}")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Generate a HomotopySolver class with the given options')
    parser.add_argument('solver_path', type=Path,
                    help='root of the solver')
    parser.add_argument('nosnoc_root', type=Path,
                    help='root of nosnoc')

    args = parser.parse_args()

    solver_path = args.solver_path
    env = Environment(loader=FileSystemLoader(str(Path(args.nosnoc_root, 'codegen_templates'))))
    
    with open(str(Path(solver_path, "solver.json"))) as f:
        solver_data = json.load(f)
        
    if not isinstance(solver_data["nlp_lbw"], list):
        solver_data["nlp_lbw"] = [solver_data["nlp_lbw"]];
    if not isinstance(solver_data["nlp_ubw"], list):
        solver_data["nlp_ubw"] = [solver_data["nlp_ubw"]];
    if not isinstance(solver_data["nlp_init_lam_w"], list):
        solver_data["nlp_init_lam_w"] = [solver_data["nlp_init_lam_w"]];
    if not isinstance(solver_data["nlp_lbg"], list):
        solver_data["nlp_lbg"] = [solver_data["nlp_lbg"]];
    if not isinstance(solver_data["nlp_ubg"], list):
        solver_data["nlp_ubg"] = [solver_data["nlp_ubg"]];
    if not isinstance(solver_data["nlp_init_lam_g"], list):
        solver_data["nlp_init_lam_g"] = [solver_data["nlp_init_lam_g"]];
    if not isinstance(solver_data["nlp_p0"], list):
        solver_data["nlp_p0"] = [solver_data["nlp_p0"]];
    if not isinstance(solver_data["nlp_x0"], list):
        solver_data["nlp_x0"] = [solver_data["nlp_x0"]];
        
    gen_HS_cpp(env, solver_data, solver_path)
    gen_HS_hpp(env, solver_data, solver_path)
    gen_cmake(env, solver_data, solver_path)

    if Path(solver_path, "problem_options.json").is_file():
        with open(str(Path(solver_path, "problem_options.json"))) as f:
            problem_options = json.load(f)
        gen_PO_hpp(env, problem_options, solver_path)

    
