def map_position(sequence: str, position: int) -> dict:
    if position < 1 or position > len(sequence):
        return {"position": position, "base": None, "valid": False}
    return {"position": position, "base": sequence[position - 1], "valid": True}


def build_position_map(sequence: str, wobble_position: int | None = 34) -> dict[int, int | None]:
    position_map: dict[int, int | None] = {}
    if wobble_position is None:
        return {position: None for position in range(1, len(sequence) + 1)}

    try:
        wobble = int(wobble_position)
    except (TypeError, ValueError):
        return {position: None for position in range(1, len(sequence) + 1)}

    offset = 34 - wobble
    for position in range(1, len(sequence) + 1):
        position_map[position] = position + offset
    return position_map
