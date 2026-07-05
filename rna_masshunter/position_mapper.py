def map_position(sequence: str, position: int) -> dict:
    if position < 1 or position > len(sequence):
        return {"position": position, "base": None, "valid": False}
    return {"position": position, "base": sequence[position - 1], "valid": True}
