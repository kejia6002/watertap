# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f"Hi, {name}")  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    print_hi("PyCharm")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
def my_decorator(func):
    def wrapper():
        print("Something is happening before the function is called.")
        func()
        print("Something is happening after the function is called.")

    return wrapper


@my_decorator
def say_hello():
    print("Hello!")


say_hello()


def repeat(n):
    def decorator(func):
        def wrapper(*args, **kwargs):
            for _ in range(n):
                func(*args, **kwargs)

        return wrapper

    return decorator


@repeat(3)
def say_hello():
    print("Hello!")


say_hello()


def flow_mol_conversion(Mw_solute, Mw_water, conc, flowrate, density):
    # flowrate in m3/s
    # conc in g/L = kg/m3
    # flow_mol in mol/s
    # density in kg/m3
    # Mw in mol/s
    flow_solute = 1e3 * conc / Mw_solute * flowrate
    flow_water = 1e3 * (flowrate * density - conc * flowrate) / Mw_water
    print("flow_solute, in molar/s", flow_solute)
    print("flow_water, in molar/s", flow_water)
    return flow_solute, flow_water


flow_mol_conversion(58.44, 18, 5, 5.2e-4, 1e3)
